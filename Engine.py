import sys
import math
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QDockWidget,
    QFormLayout, QLineEdit, QPushButton
)
from PyQt5.QtGui import QPainter, QPen, QColor, QFont, QPolygonF
from PyQt5.QtCore import Qt, QTimer, QPointF


SAFE_GLOBALS = {
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "sqrt": math.sqrt,
    "pi": math.pi,
    "abs": abs,
    "log": math.log,
    "exp": math.exp,
    "pow": math.pow
}


class GeometryCalculator:
    def __init__(self, eq_x, eq_y, eq_z,
                 u_min, u_max, u_steps,
                 v_min, v_max, v_steps,
                 a_values, b_values):
        self.eq_x = eq_x
        self.eq_y = eq_y
        self.eq_z = eq_z
        self.u_min = u_min
        self.u_max = u_max
        self.u_steps = u_steps
        self.v_min = v_min
        self.v_max = v_max
        self.v_steps = v_steps
        self.a_values = a_values
        self.b_values = b_values
        self.compile_equations()

    def compile_equations(self):
        try:
            self.func_x = lambda u, v, a, b: eval(
                self.eq_x, SAFE_GLOBALS,
                {"u": u, "v": v, "a": a, "b": b}
            )
            self.func_y = lambda u, v, a, b: eval(
                self.eq_y, SAFE_GLOBALS,
                {"u": u, "v": v, "a": a, "b": b}
            )
            self.func_z = lambda u, v, a, b: eval(
                self.eq_z, SAFE_GLOBALS,
                {"u": u, "v": v, "a": a, "b": b}
            )
        except Exception as e:
            print("EQUATIONS COMPILATION ERROR", e)

    def make_range(self, start, stop, step):
        values = []
        v = start
        while v <= stop + 1e-9:
            values.append(round(v, 10))
            v += step
        return values

    def compute_points_and_normals(self, a, b):
        points = []
        try:
            du = (self.u_max - self.u_min) / (self.u_steps - 1)
            dv = (self.v_max - self.v_min) / (self.v_steps - 1)
        except Exception as e:
            print("DU OR DV ERROR:", e)
            return points, []
        for i in range(self.u_steps):
            row = []
            u = self.u_min + i * du
            for j in range(self.v_steps):
                v = self.v_min + j * dv
                try:
                    x = self.func_x(u, v, a, b)
                    y = self.func_y(u, v, a, b)
                    z = self.func_z(u, v, a, b)
                except Exception as e:
                    print(f"POINT CALCULATION ERROR (u={u}, v={v}, a={a}, b={b}):", e)
                    x = y = z = 0
                row.append((x, y, z))
            points.append(row)

        normals = [[(0, 0, 0) for _ in range(self.v_steps)] for _ in range(self.u_steps)]
        for i in range(self.u_steps):
            for j in range(self.v_steps):
                try:
                    if i < self.u_steps - 1 and j < self.v_steps - 1:
                        p1 = points[i][j]
                        p2 = points[i+1][j]
                        p3 = points[i][j+1]
                        normals[i][j] = self.surface_normal(p1, p2, p3)
                    else:
                        normals[i][j] = (0, 0, 1)
                except Exception as e:
                    print(f"NORMAL CALCULATION ERROR for point ({i},{j}):", e)
                    normals[i][j] = (0, 0, 1)
        return points, normals

    def surface_normal(self, p1, p2, p3):
        try:
            ux = p2[0] - p1[0]
            uy = p2[1] - p1[1]
            uz = p2[2] - p1[2]
            vx = p3[0] - p1[0]
            vy = p3[1] - p1[1]
            vz = p3[2] - p1[2]
            nx = uy * vz - uz * vy
            ny = uz * vx - ux * vz
            nz = ux * vy - uy * vx
            length = math.sqrt(nx*nx + ny*ny + nz*nz)
            if length == 0:
                return (0, 0, 0)
            return (nx/length, ny/length, nz/length)
        except Exception as e:
            print("SURFACE NORMAL CALCULATION ERROR:", e)
            return (0, 0, 0)

    def compute_triangles(self):
        triangles = []
        try:
            for a in self.a_values:
                for b in self.b_values:
                    points, normals = self.compute_points_and_normals(a, b)
                    for i in range(self.u_steps - 1):
                        for j in range(self.v_steps - 1):
                            try:
                                tri1 = {
                                    'vertices': [points[i][j], points[i+1][j], points[i][j+1]],
                                    'normals': [normals[i][j], normals[i+1][j], normals[i][j+1]]
                                }
                                triangles.append(tri1)
                            except Exception as e:
                                print(f"FIRST TRIANGLE CREATION ERROR for cell ({i},{j}):", e)
                            try:
                                tri2 = {
                                    'vertices': [points[i+1][j], points[i+1][j+1], points[i][j+1]],
                                    'normals': [normals[i+1][j], normals[i+1][j+1], normals[i][j+1]]
                                }
                                triangles.append(tri2)
                            except Exception as e:
                                print(f"SECOND TRIANGLE CREATION ERROR for cell ({i},{j}):", e)
        except Exception as e:
            print("TRIANGLE COMPUTATION ERROR:", e)
        return triangles


class GraphWidget(QWidget):
    def __init__(self):
        super().__init__()

        self.eq_x = "(a + b*cos(v))*cos(u)"
        self.eq_y = "(a + b*cos(v))*sin(u)"
        self.eq_z = "b*sin(v) + a*u"

        self.a_min = 1.5
        self.a_max = 1.5
        self.a_step = 1.0
        self.b_min = 0.5
        self.b_max = 0.5
        self.b_step = 1.0

        self.u_min = 0.0
        self.u_max = 4 * math.pi
        self.u_steps = 30    

        self.v_min = 0.0
        self.v_max = 2 * math.pi
        self.v_steps = 20   

        a_values = self.make_range(self.a_min, self.a_max, self.a_step)
        b_values = self.make_range(self.b_min, self.b_max, self.b_step)
        self.geometry = GeometryCalculator(
            self.eq_x, self.eq_y, self.eq_z,
            self.u_min, self.u_max, self.u_steps,
            self.v_min, self.v_max, self.v_steps,
            a_values, b_values
        )


        try:
            self.triangles = self.geometry.compute_triangles()
        except Exception as e:
            print("TRIANGLE CACHING ERROR:", e)
            self.triangles = []

        self.rotation_x = -math.pi / 60
        self.rotation_z = 0.0
        self.zoom = 2.0     
        self.translation_x = 0
        self.translation_y = 0
        self.camera_distance = 100
        self.light_direction = self.normalize((0, -1, 1))
        self.last_mouse_pos = None
        self.top_view_toggle = False 

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update)
        self.timer.start(30)
        self.setFocusPolicy(Qt.StrongFocus)

    def make_range(self, start, stop, step):
        values = []
        v = start
        while v <= stop + 1e-9:
            values.append(round(v, 10))
            v += step
        return values

    def normalize(self, vec):
        try:
            length = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
            if length == 0:
                return (0, 0, 0)
            return (vec[0]/length, vec[1]/length, vec[2]/length)
        except Exception as e:
            print("NORMALIZATION ERROR:", e)
            return (0, 0, 0)

    def project_point(self, point):
        try:
            x, y, z = point
            cosx = math.cos(self.rotation_x)
            sinx = math.sin(self.rotation_x)
            y2 = y * cosx - z * sinx
            z2 = y * sinx + z * cosx
            x2 = x
       
            cosz = math.cos(self.rotation_z)
            sinz = math.sin(self.rotation_z)
            x3 = x2 * cosz - y2 * sinz
            y3 = x2 * sinz + y2 * cosz
            z3 = z2

            d = self.camera_distance
            factor = d / (d + y3) if (d + y3) != 0 else 1
            screen_x = x3 * factor * self.zoom
            screen_y = z3 * factor * self.zoom
            center_x = self.width() / 2 + self.translation_x
            center_y = self.height() / 2 + self.translation_y
            return (center_x + screen_x, center_y - screen_y)
        except Exception as e:
            print("PROJECT POINT ERROR for point", point, ":", e)
            return (0, 0)

    def dot(self, a, b):
        try:
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
        except Exception as e:
            print("DOT PRODUCT ERROR for vectors", a, b, ":", e)
            return 0

    def paintEvent(self, event):
        try:
            painter = QPainter(self)
            painter.setRenderHint(QPainter.Antialiasing)
            painter.fillRect(self.rect(), QColor(255, 255, 255))
        except Exception as e:
            print("QPAINTER INITIALIZATION ERROR:", e)
            return

        try:
            self.draw_grid(painter)
        except Exception as e:
            print("GRID DRAWING ERROR:", e)

        try:
           
            sorted_triangles = sorted(
                self.triangles,
                key=lambda tri: sum(self.get_depth(pt) for pt in tri['vertices']) / 3,
                reverse=True
            )

            for tri in sorted_triangles:
                try:
                    vertices = tri['vertices']
                    normals = tri['normals']
                    projected = [self.project_point(pt) for pt in vertices]
                    rotated_normals = [self.rotate_vector(n) for n in normals]
                    brightness = sum(max(0, self.dot(n, self.light_direction)) for n in rotated_normals) / 3
                    color_intensity = int(50 + 205 * brightness)
                    fill_color = QColor(100, color_intensity, color_intensity)
                    poly = QPolygonF([QPointF(x, y) for (x, y) in projected])
                    painter.setBrush(fill_color)
                    painter.drawPolygon(poly)
                    painter.setBrush(Qt.NoBrush)
                    pen = QPen(QColor(0, 0, 0))
                    pen.setWidth(1)
                    painter.setPen(pen)
                    painter.drawLine(QPointF(*projected[0]), QPointF(*projected[1]))
                    painter.drawLine(QPointF(*projected[1]), QPointF(*projected[2]))
                    painter.drawLine(QPointF(*projected[2]), QPointF(*projected[0]))
                except Exception as e:
                    print("TRIANGLE DRAWING ERROR:", e)
        except Exception as e:
            print("TRIANGLES RENDERING ERROR:", e)

        try:
            self.draw_axes(painter)
        except Exception as e:
            print("AXES DRAWING ERROR:", e)

    def draw_grid(self, painter):
        try:
            grid_spacing = 50
            pen = QPen(QColor(200, 200, 200))
            pen.setStyle(Qt.DotLine)
            painter.setPen(pen)
            w = self.width()
            h = self.height()
            for x in range(0, w, grid_spacing):
                painter.drawLine(x, 0, x, h)
            for y in range(0, h, grid_spacing):
                painter.drawLine(0, y, w, y)
        except Exception as e:
            print("GRID DRAWING ERROR:", e)

    def draw_axes(self, painter):
        try:
            center_x = self.width() / 2 + self.translation_x
            center_y = self.height() / 2 + self.translation_y
            origin = (center_x, center_y)
            p_x = self.project_point((100, 0, 0))
            p_y = self.project_point((0, 100, 0))
            p_z = self.project_point((0, 0, 100))
            pen = QPen()
            pen.setWidth(2)
            pen.setColor(QColor(255, 0, 0))
            painter.setPen(pen)
            painter.drawLine(int(origin[0]), int(origin[1]), int(p_x[0]), int(p_x[1]))
            pen.setColor(QColor(0, 255, 0))
            painter.setPen(pen)
            painter.drawLine(int(origin[0]), int(origin[1]), int(p_y[0]), int(p_y[1]))
            pen.setColor(QColor(0, 0, 255))
            painter.setPen(pen)
            painter.drawLine(int(origin[0]), int(origin[1]), int(p_z[0]), int(p_z[1]))
            painter.setPen(QPen(Qt.black))
            painter.setFont(QFont("Arial", 10))
            painter.drawText(int(p_x[0]), int(p_x[1]), "X")
            painter.drawText(int(p_y[0]), int(p_y[1]), "Y")
            painter.drawText(int(p_z[0]), int(p_z[1]), "Z")
        except Exception as e:
            print("AXES DRAWING ERROR:", e)

    def keyPressEvent(self, event):
        try:
            if event.key() == Qt.Key_Up:
                self.top_view_toggle = not self.top_view_toggle
                if self.top_view_toggle:
                    self._orig_rotation_x = self.rotation_x
                    self.rotation_x = -math.pi / 2
                else:
                    if hasattr(self, "_orig_rotation_x"):
                        self.rotation_x = self._orig_rotation_x
                self.update()
            else:
                super().keyPressEvent(event)
        except Exception as e:
            print("KEY PRESS EVENT ERROR:", e)

    def mousePressEvent(self, event):
        try:
            if event.button() == Qt.LeftButton:
                self.last_mouse_pos = event.pos()
        except Exception as e:
            print("MOUSE PRESS EVENT ERROR:", e)

    def mouseMoveEvent(self, event):
        try:
            if event.buttons() & Qt.LeftButton and self.last_mouse_pos is not None:
                delta = event.pos() - self.last_mouse_pos
                self.rotation_z += delta.x() * 0.01
                self.last_mouse_pos = event.pos()
                self.update()
        except Exception as e:
            print("MOUSE MOVE EVENT ERROR:", e)

    def mouseReleaseEvent(self, event):
        try:
            self.last_mouse_pos = None
        except Exception as e:
            print("MOUSE RELEASE EVENT ERROR:", e)

    def wheelEvent(self, event):
        try:
            delta = event.angleDelta().y()
            if delta > 0:
                self.zoom *= 1.1
            else:
                self.zoom /= 1.1
            self.update()
        except Exception as e:
            print("WHEEL EVENT ERROR:", e)

    def update_parameters(self, params):
        try:
            self.eq_x = params.get("eq_x", self.eq_x)
            self.eq_y = params.get("eq_y", self.eq_y)
            self.eq_z = params.get("eq_z", self.eq_z)
            self.a_min = float(params.get("a_min", self.a_min))
            self.a_max = float(params.get("a_max", self.a_max))
            self.a_step = float(params.get("a_step", self.a_step))
            self.b_min = float(params.get("b_min", self.b_min))
            self.b_max = float(params.get("b_max", self.b_max))
            self.b_step = float(params.get("b_step", self.b_step))
            self.u_min = float(params.get("u_min", self.u_min))
            self.u_max = float(params.get("u_max", self.u_max))
            self.u_steps = int(params.get("u_steps", self.u_steps))
            self.v_min = float(params.get("v_min", self.v_min))
            self.v_max = float(params.get("v_max", self.v_max))
            self.v_steps = int(params.get("v_steps", self.v_steps))

            a_values = self.make_range(self.a_min, self.a_max, self.a_step)
            b_values = self.make_range(self.b_min, self.b_max, self.b_step)
            self.geometry = GeometryCalculator(
                self.eq_x, self.eq_y, self.eq_z,
                self.u_min, self.u_max, self.u_steps,
                self.v_min, self.v_max, self.v_steps,
                a_values, b_values
            )
            self.triangles = self.geometry.compute_triangles()
            self.update()
        except Exception as e:
            print("PARAMETER UPDATE ERROR:", e)

    def rotate_vector(self, vec):
        x, y, z = vec
        cosx = math.cos(self.rotation_x)
        sinx = math.sin(self.rotation_x)
        y2 = y * cosx - z * sinx
        z2 = y * sinx + z * cosx
        x2 = x
    
        cosz = math.cos(self.rotation_z)
        sinz = math.sin(self.rotation_z)
        x3 = x2 * cosz - y2 * sinz
        y3 = x2 * sinz + y2 * cosz
        return (x3, y3, z2)
    
    def get_depth(self, point):
        x, y, z = point
        cosx = math.cos(self.rotation_x)
        sinx = math.sin(self.rotation_x)
        y2 = y * cosx - z * sinx
        x2 = x
        cosz = math.cos(self.rotation_z)
        sinz = math.sin(self.rotation_z)
        y3 = x2 * sinz + y2 * cosz
        return y3


class ControlPanel(QWidget):
    def __init__(self, graph_widget):
        super().__init__()
        self.graph_widget = graph_widget
        self.init_ui()

    def init_ui(self):
        try:
            layout = QFormLayout()
            self.le_eq_x = QLineEdit(self.graph_widget.eq_x)
            self.le_eq_y = QLineEdit(self.graph_widget.eq_y)
            self.le_eq_z = QLineEdit(self.graph_widget.eq_z)
            layout.addRow("x(u,v) =", self.le_eq_x)
            layout.addRow("y(u,v) =", self.le_eq_y)
            layout.addRow("z(u,v) =", self.le_eq_z)
            self.le_a_min = QLineEdit(str(self.graph_widget.a_min))
            self.le_a_max = QLineEdit(str(self.graph_widget.a_max))
            self.le_a_step = QLineEdit(str(self.graph_widget.a_step))
            layout.addRow("a min =", self.le_a_min)
            layout.addRow("a max =", self.le_a_max)
            layout.addRow("a step =", self.le_a_step)
            self.le_b_min = QLineEdit(str(self.graph_widget.b_min))
            self.le_b_max = QLineEdit(str(self.graph_widget.b_max))
            self.le_b_step = QLineEdit(str(self.graph_widget.b_step))
            layout.addRow("b min =", self.le_b_min)
            layout.addRow("b max =", self.le_b_max)
            layout.addRow("b step =", self.le_b_step)
            self.le_u_min = QLineEdit(str(self.graph_widget.u_min))
            self.le_u_max = QLineEdit(str(self.graph_widget.u_max))
            self.le_u_steps = QLineEdit(str(self.graph_widget.u_steps))
            layout.addRow("u min =", self.le_u_min)
            layout.addRow("u max =", self.le_u_max)
            layout.addRow("u steps =", self.le_u_steps)
            self.le_v_min = QLineEdit(str(self.graph_widget.v_min))
            self.le_v_max = QLineEdit(str(self.graph_widget.v_max))
            self.le_v_steps = QLineEdit(str(self.graph_widget.v_steps))
            layout.addRow("v min =", self.le_v_min)
            layout.addRow("v max =", self.le_v_max)
            layout.addRow("v steps =", self.le_v_steps)
            self.btn_update = QPushButton("Обновить")
            self.btn_update.clicked.connect(self.on_update)
            layout.addRow(self.btn_update)
            self.setLayout(layout)
        except Exception as e:
            print("CONTROL PANEL CREATION ERROR:", e)

    def on_update(self):
        try:
            params = {
                "eq_x": self.le_eq_x.text(),
                "eq_y": self.le_eq_y.text(),
                "eq_z": self.le_eq_z.text(),
                "a_min": self.le_a_min.text(),
                "a_max": self.le_a_max.text(),
                "a_step": self.le_a_step.text(),
                "b_min": self.le_b_min.text(),
                "b_max": self.le_b_max.text(),
                "b_step": self.le_b_step.text(),
                "u_min": self.le_u_min.text(),
                "u_max": self.le_u_max.text(),
                "u_steps": self.le_u_steps.text(),
                "v_min": self.le_v_min.text(),
                "v_max": self.le_v_max.text(),
                "v_steps": self.le_v_steps.text()
            }
            self.graph_widget.update_parameters(params)
        except Exception as e:
            print("CONTROL PANEL UPDATE ERROR:", e)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        try:
            self.setWindowTitle("Graphic Interface for Functions")
            self.graph_widget = GraphWidget()
            self.setCentralWidget(self.graph_widget)
            self.control_panel = ControlPanel(self.graph_widget)
            self.dock = QDockWidget("Параметры", self)
            self.dock.setWidget(self.control_panel)
            self.dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
            self.addDockWidget(Qt.LeftDockWidgetArea, self.dock)
            self.dock.setFixedWidth(300)
            self.resize(1000, 700)
        except Exception as e:
            print("MAIN WINDOW INITIALIZATION ERROR:", e)

if __name__ == '__main__':
    try:
        app = QApplication(sys.argv)
        window = MainWindow()
        window.show()
        sys.exit(app.exec_())
    except Exception as e:
        print("MAIN APPLICATION LOOP ERROR:", e)
