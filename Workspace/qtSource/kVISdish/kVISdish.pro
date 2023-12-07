######################################################################
# Automatically generated by qmake (3.1) Wed Dec 6 09:47:08 2023
######################################################################

TEMPLATE = app
TARGET = kVISdish
INCLUDEPATH += .

QT += widgets

# You can make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# Please consult the documentation of the deprecated API in order to know
# how to port your code away from it.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input
HEADERS += kVISdish.h \

FORMS += kVISdish.ui
SOURCES += kVISdish.cpp \
           main.cpp \


CONFIG += link_pkgconfig
PKGCONFIG += opencv4

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../Temp/Vimba_6_0/VimbaCPP/DynamicLib/x86_64bit/release/ -lVimbaCPP
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../Temp/Vimba_6_0/VimbaCPP/DynamicLib/x86_64bit/debug/ -lVimbaCPP
else:unix: LIBS += -L$$PWD/../../../Temp/Vimba_6_0/VimbaCPP/DynamicLib/x86_64bit/ -lVimbaCPP

INCLUDEPATH += $$PWD/../../../Temp/Vimba_6_0/VimbaCPP/Include
DEPENDPATH += $$PWD/../../../Temp/Vimba_6_0/VimbaCPP/Include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../Temp/Vimba_6_0/VimbaC/DynamicLib/x86_64bit/release/ -lVimbaC
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../Temp/Vimba_6_0/VimbaC/DynamicLib/x86_64bit/debug/ -lVimbaC
else:unix: LIBS += -L$$PWD/../../../Temp/Vimba_6_0/VimbaC/DynamicLib/x86_64bit/ -lVimbaC

INCLUDEPATH += $$PWD/../../../Temp/Vimba_6_0/VimbaC/Include
DEPENDPATH += $$PWD/../../../Temp/Vimba_6_0/VimbaC/Include