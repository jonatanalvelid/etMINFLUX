from PyQt5 import QtWidgets, QtCore, QtGui

class PushButton(QtWidgets.QPushButton):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.setStyleSheet('background-color: rgb(50,50,50); color: rgb(217,83,0);')

class FieldLabel(QtWidgets.QLabel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.setStyleSheet('color: rgb(217,83,0);')

class TitleLabel(QtWidgets.QLabel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold, pointSize=12))
        self.setStyleSheet('background-color: rgb(30,30,30); color: rgb(217,83,0);')

class ComboBox(QtWidgets.QComboBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setStyleSheet('background-color: rgb(50,50,50); color: rgb(170,170,170); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')

class CheckBox(QtWidgets.QCheckBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setStyleSheet('background-color: transparent; color: rgb(217,83,0);')

class LineEdit(QtWidgets.QLineEdit):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stylesheet = 'background-color: rgb(50,50,50); color: rgb(170,170,170); border: 1px solid rgb(100,100,100);'
        self.setStyleSheet(self.stylesheet)

    def setEditable(self, editable):
        if editable:
            self.stylesheet = 'background-color: rgb(50,50,50); color: rgb(170,170,170); border: 1px solid rgb(100,100,100);'
            self.setReadOnly(False)
        else:
            self.stylesheet = 'background-color: rgb(50,50,50); color: rgb(100,100,100); border: 1px solid rgb(100,100,100);'
            self.setReadOnly(True)
        self.setStyleSheet(self.stylesheet)
        

