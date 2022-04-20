import sys

from PyQt6.QtWidgets import QMainWindow, QApplication
from UI.MainWindow import Ui_MainWindow


class MainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.setWindowTitle("Mr. Chem Gui")

        self.button_quit.clicked.connect(self.close_button_clicked)

    def close_button_clicked(self):
        self.close()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()