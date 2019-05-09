#-------------------------------------------------
#
# Project created by QtCreator 2019-04-17T10:42:45
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = tempFEM
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        main.cpp \
        widget.cpp \
    temp2dfemcore.cpp \
    qcustomplot/qcustomplot.cpp \
#    mainwindow.cpp
    2DInterpolation/interpolation.cpp \
    metis-5.1.0/programs/mpmetis.c \
    metis-5.1.0/programs/cmdline_gpmetis.c \
    metis-5.1.0/programs/io.c

HEADERS += \
        widget.h \
    temp2dfemcore.h \
    datatype.h \
    qcustomplot/qcustomplot.h \
    metis-5.1.0/programs/mpmetis.h
#    mainwindow.h

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../armadillo/examples/lib_win64/ -lblas_win64_MT
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../armadillo/examples/lib_win64/ -lblas_win64_MT
else:unix: LIBS += -L$$PWD/../../../armadillo/examples/lib_win64/ -lblas_win64_MT

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../armadillo/examples/lib_win64/ -llapack_win64_MT
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../armadillo/examples/lib_win64/ -llapack_win64_MT
else:unix: LIBS += -L$$PWD/../../../armadillo/examples/lib_win64/ -llapack_win64_MT

INCLUDEPATH += $$PWD/../../../armadillo/include

DEPENDPATH += $$PWD/../../../armadillo/include

CONFIG +=console

#FORMS += \
#    mainwindow.ui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

INCLUDEPATH +=     metis-5.1.0/GKlib\
    metis-5.1.0/include \
    metis-5.1.0/programs \
    SuperLU_5.2.1/SRC \
    SuperLU_5.2.1/CBLAS \


LIBS += D:\tempFEM\tempFEM0\tempFEM\metis.lib\
        D:\DDTLM\SuperLU_5.2.1\SuperLU\x64\Release\SuperLU.lib \
        D:\DDTLM\SuperLU_5.2.1\SuperLU\x64\Release\CBLAS.lib \
