######################################################################
# Automatically generated by qmake (3.1) Tue Dec 5 08:14:06 2023
######################################################################

TEMPLATE = app
TARGET = SimulDish
INCLUDEPATH += .

QT += widgets

# You can make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# Please consult the documentation of the deprecated API in order to know
# how to port your code away from it.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input
HEADERS += SimulDish.h
FORMS += SimulDish.ui
SOURCES += main.cpp SimulDish.cpp

CONFIG += link_pkgconfig
PKGCONFIG += opencv4 qwt
