######################################################################
# Automatically generated by qmake (3.1) Tue Jun 27 09:38:56 2023
######################################################################

TEMPLATE = app
TARGET = kFluxMapper
INCLUDEPATH += .

QT += widgets

# You can make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# Please consult the documentation of the deprecated API in order to know
# how to port your code away from it.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input
HEADERS += kfluxmapper.h
FORMS += kfluxmapper.ui
SOURCES += kfluxmapper.cpp main.cpp

CONFIG += link_pkgconfig
PKGCONFIG += opencv4 glib-2.0 aravis-0.6
