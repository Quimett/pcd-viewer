import numpy as np
import os
import time
import random
import sys

try:
    import pyvista as pv
except ImportError:
    pv = None

class LectorPCD:
    """
    Lee un archivo PCD ASCII y extrae las coordenadas x, y, z.
    """
    @staticmethod
    def leer_pcd(ruta_archivo):
        puntos = []
        with open(ruta_archivo, 'r') as f:
            for linea in f:
                if linea.strip().lower() == 'data ascii':
                    break
            for linea in f:
                partes = linea.strip().split()
                if len(partes) < 3:
                    continue
                x, y, z = map(float, partes[:3])
                puntos.append((x, y, z))
        return np.array(puntos)

class RejillaOcupacion:
    """
    Cuadrícula de ocupación uniforme que almacena recuentos por celda cúbica.
    """
    def __init__(self, puntos, tam_celda):
        self.puntos = puntos
        self.tam_celda = tam_celda
        self.x_minimo, self.y_minimo, self.z_minimo = np.min(puntos, axis=0)
        self.x_maximo, self.y_maximo, self.z_maximo = np.max(puntos, axis=0)
        nx = int(np.floor((self.x_maximo - self.x_minimo) / tam_celda)) + 1
        ny = int(np.floor((self.y_maximo - self.y_minimo) / tam_celda)) + 1
        nz = int(np.floor((self.z_maximo - self.z_minimo) / tam_celda)) + 1
        self.total_celdas = nx * ny * nz
        self.rejilla = {}
        self._poblar()

    def _poblar(self):
        for x, y, z in self.puntos:
            ix = int((x - self.x_minimo) // self.tam_celda)
            iy = int((y - self.y_minimo) // self.tam_celda)
            iz = int((z - self.z_minimo) // self.tam_celda)
            clave = (ix, iy, iz)
            self.rejilla[clave] = self.rejilla.get(clave, 0) + 1

    def estadisticas_celdas(self):
        ocupadas = len(self.rejilla)
        vacias = self.total_celdas - ocupadas
        promedio_puntos = np.mean(list(self.rejilla.values())) if ocupadas else 0
        return {'total_celdas': self.total_celdas, 'ocupadas': ocupadas, 'vacias': vacias, 'promedio_puntos': promedio_puntos}

class NodoOctree:
    """
    Octree adaptable: subdivide nodos con > maximo_puntos hasta que el tamaño <= tam_minimo.
    """
    def __init__(self, puntos, limites, tam_minimo, maximo_puntos):
        self.puntos = puntos
        self.limites = limites
        self.tam_minimo = tam_minimo
        self.maximo_puntos = maximo_puntos
        self.hijos = []
        self._subdividir_o_hoja()

    def _subdividir_o_hoja(self):
        tam = np.subtract(self.limites[1], self.limites[0])
        if len(self.puntos) > self.maximo_puntos and np.max(tam) > self.tam_minimo:
            self._subdividir()

    def _subdividir(self):
        (x0, y0, z0), (x1, y1, z1) = self.limites
        mx, my, mz = (x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2
        for dx in [(x0, mx), (mx, x1)]:
            for dy in [(y0, my), (my, y1)]:
                for dz in [(z0, mz), (mz, z1)]:
                    sub_limites = ((dx[0], dy[0], dz[0]), (dx[1], dy[1], dz[1]))
                    mascara = (
                        (self.puntos[:,0] >= sub_limites[0][0]) & (self.puntos[:,0] < sub_limites[1][0]) &
                        (self.puntos[:,1] >= sub_limites[0][1]) & (self.puntos[:,1] < sub_limites[1][1]) &
                        (self.puntos[:,2] >= sub_limites[0][2]) & (self.puntos[:,2] < sub_limites[1][2])
                    )
                    puntos_hijo = self.puntos[mascara]
                    self.hijos.append(NodoOctree(puntos_hijo, sub_limites, self.tam_minimo, self.maximo_puntos))

    def recopilar_estadisticas(self):
        def recursivo(nodo):
            estad = {'nodos': 1, 'hojas': 0, 'internos': 0, 'hojas_ocupadas': 0,
                     'hojas_vacias': 0, 'suma_puntos': 0}
            if not nodo.hijos:
                estad['hojas'] = 1
                cnt = len(nodo.puntos)
                if cnt > 0:
                    estad['hojas_ocupadas'] = 1
                    estad['suma_puntos'] = cnt
                else:
                    estad['hojas_vacias'] = 1
            else:
                estad['internos'] = 1
                for h in nodo.hijos:
                    sub = recursivo(h)
                    for k in estad:
                        estad[k] += sub[k]
            return estad
        raw = recursivo(self)
        prom_hojas = raw['suma_puntos'] / raw['hojas_ocupadas'] if raw['hojas_ocupadas'] else 0
        return {'total_nodos': raw['nodos'], 'hojas': raw['hojas'], 'internos': raw['internos'],
                'hojas_ocupadas': raw['hojas_ocupadas'], 'hojas_vacias': raw['hojas_vacias'], 'promedio_puntos': prom_hojas}

    def obtener_nodos_hoja(self):
        if not self.hijos:
            return [self]
        hojas = []
        for h in self.hijos:
            hojas.extend(h.obtener_nodos_hoja())
        return hojas

PARAMETROS = [
    {'tam_celda': 1.0, 'tam_minimo': 1.0, 'maximo_puntos': 100},
    {'tam_celda': 0.5, 'tam_minimo': 0.5, 'maximo_puntos': 100},
    {'tam_celda': 3.0, 'tam_minimo': 1.0, 'maximo_puntos': 300},
]

def comparar_archivo(ruta):
    pts = LectorPCD.leer_pcd(ruta)
    pts = np.unique(pts, axis=0)
    limites = ((pts[:,0].min(), pts[:,1].min(), pts[:,2].min()),
               (pts[:,0].max(), pts[:,1].max(), pts[:,2].max()))
    resultados = []
    for p in PARAMETROS:
        g = RejillaOcupacion(pts, p['tam_celda']).estadisticas_celdas()
        o = NodoOctree(pts, limites, p['tam_minimo'], p['maximo_puntos']).recopilar_estadisticas()
        resultados.append((p, g, o))
    return resultados

def formatear_comparacion(nombre, datos):
    columnas = []
    for p, g, o in datos:
        columnas.append([
            f"Analizando: {nombre}",
            "Análisis Comparativo",
            f"Archivo: {nombre}",
            "Rejilla:",
            f"  tam_celda: {p['tam_celda']}",
            f"  total_celdas: {g['total_celdas']}",
            f"  ocupadas: {g['ocupadas']}",
            f"  vacias: {g['vacias']}",
            f"  promedio_puntos: {g['promedio_puntos']:.2f}",
            "Octree:",
            f"  tam_minimo: {p['tam_minimo']}, maximo_puntos: {p['maximo_puntos']}",
            f"  total_nodos: {o['total_nodos']}",
            f"  hojas: {o['hojas']}",
            f"  hojas_ocupadas: {o['hojas_ocupadas']}",
            f"  hojas_vacias: {o['hojas_vacias']}",
            f"  promedio_puntos: {o['promedio_puntos']:.2f}"
        ])
    L = max(len(c) for c in columnas)
    W = max(max(len(x) for x in c) for c in columnas) + 2
    separador = ' | '
    borde = '-' * (len(columnas)*W + (len(columnas)-1)*len(separador))
    print(borde)
    for i in range(L):
        linea = ''.join(columnas[j][i].ljust(W) + (separador if j < len(columnas)-1 else '')
                        for j in range(len(columnas)))
        print(linea)
    print(borde)

def visualizar_archivo(ruta):
    if pv is None:
        print("Error: PyVista no está instalado. Instale con 'pip install pyvista'.")
        sys.exit(1)
    pts = LectorPCD.leer_pcd(ruta)
    print(f"Cargados {len(pts)} puntos.")
    limites = ((pts[:,0].min(), pts[:,1].min(), pts[:,2].min()),
               (pts[:,0].max(), pts[:,1].max(), pts[:,2].max()))
    dims = np.ptp(pts, axis=0)
    max_dim = max(dims)
    tam_minimo = max_dim / 50.0
    total_puntos = pts.shape[0]
    maximo_puntos = max(int(total_puntos / 1000), 10)
    max_hojas = min(int(total_puntos / 10), 10000)
    print(f"Parámetros auto -> tam_minimo: {tam_minimo:.3f}, maximo_puntos: {maximo_puntos}, max_hojas: {max_hojas}")
    print("Construyendo Octree...")
    t0 = time.time()
    octree = NodoOctree(pts, limites, tam_minimo, maximo_puntos)
    print(f"Construido en {time.time()-t0:.2f}s")
    hojas = octree.obtener_nodos_hoja()
    hojas_ocupadas = [n for n in hojas if len(n.puntos) > 0]
    T = len(hojas_ocupadas)
    print(f"Hojas ocupadas: {T}")
    if max_hojas > 0 and T > max_hojas:
        print(f"Seleccionando {max_hojas} hojas representativas...")
        centroides = np.array([
            ((b[0][0]+b[1][0])/2, (b[0][1]+b[1][1])/2, (b[0][2]+b[1][2])/2)
            for b in [n.limites for n in hojas_ocupadas]
        ])
        paso = T / max_hojas
        seleccionadas = [hojas_ocupadas[int(i*paso)] for i in range(max_hojas)]
    else:
        seleccionadas = hojas_ocupadas
        print(f"Dibujando todas las {T} hojas.")
    pl = pv.Plotter()
    pl.add_mesh(pv.PolyData(pts), color='black', point_size=2)
    for nodo in seleccionadas:
        (x0, y0, z0), (x1, y1, z1) = nodo.limites
        centro = ((x0+x1)/2, (y0+y1)/2, (z0+z1)/2)
        longitudes = (x1-x0, y1-y0, z1-z0)
        cubo = pv.Cube(center=centro, x_length=longitudes[0], y_length=longitudes[1], z_length=longitudes[2])
        pl.add_mesh(cubo, style='wireframe', color='green', opacity=0.4)
    pl.add_axes()
    pl.show()

def seleccionar_archivo_unico():
    archivos = [f for f in os.listdir('.') if f.lower().endswith('.pcd')]
    for i, f in enumerate(archivos, 1):
        print(f"  {i}) {f}")
    opcion = input(f"Selecciona un archivo [1-{len(archivos)}]: ")
    try:
        idx = int(opcion) - 1
        if 0 <= idx < len(archivos):
            return archivos[idx]
    except ValueError:
        pass
    print("Opción inválida.")
    sys.exit(1)

if __name__ == '__main__':
    print("=== Ejercicio 1 ===")
    print("1) Comparar Datos (Rejilla vs Octree)")
    print("2) Visualizar Octree 3D")
    modo = input("Introduce opción [1-2]: ")
    print()
    if modo == '1':
        print("-> Modo Comparar Datos")
        archivos = [f for f in os.listdir('src') if f.lower().endswith('.pcd')]
        for idx, fn in enumerate(archivos, 1):
            print(f"  {idx}) {fn}")
        print("  0) Todos los archivos")
        opcion = input(f"Selecciona archivo [0-{len(archivos)}]: ")
        if opcion == '0':
            seleccionados = archivos
        else:
            try:
                idx = int(opcion) - 1
                if idx < 0 or idx >= len(archivos):
                    raise ValueError
                seleccionados = [archivos[idx]]
            except ValueError:
                print("Opción inválida. Saliendo.")
                sys.exit(1)
        for f in seleccionados:
            ruta = os.path.join('src', f)
            print()
            formatear_comparacion(f, comparar_archivo(ruta))

    elif modo == '2':
        print("-> Modo Visualización 3D")
        archivos = [f for f in os.listdir('src') if f.lower().endswith('.pcd')]
        for i, f in enumerate(archivos, 1):
            print(f"  {i}) {f}")
        opcion = input(f"Selecciona un archivo [1-{len(archivos)}]: ")
        try:
            idx = int(opcion) - 1
            if 0 <= idx < len(archivos):
                ruta = os.path.join('src', archivos[idx])
                print()
                visualizar_archivo(ruta)
            else:
                raise ValueError
        except ValueError:
            print("Opción inválida.")
            sys.exit(1)
    else:
        print("Opción inválida. Saliendo.")