QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += main.cpp \
    ffasettings.cpp \
    ffa.cpp \
    ffa_tools.cpp \
    ffa_stages.cpp \
    testsignal.cpp \
    myqvector.cpp \
    ffago.cpp \
    ffatools.cpp \
    ffastages.cpp

HEADERS += \
    ffasettings.h \
    ffa.h \
    ffa_tools.h \
    ffa_stages.h \
    testsignal.h \
    myqvector.h \
    ffago.h \
    ffatools.h \
    ffastages.h \
    mymath.h
