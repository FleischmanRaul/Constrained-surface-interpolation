QT += core gui widgets opengl

unix {
    # for GLEW installed into /usr/lib/libGLEW.so or /usr/lib/glew.lib
    # if libraries GLU or OpenGL are required uncomment the commented libraries

     LIBS += -lGLEW -lglut -lGLU
}

FORMS += \
    GUI/MainWindow.ui \
    GUI/SideWidget.ui

HEADERS += \
    GUI/GLWidget.h \
    GUI/MainWindow.h \
    GUI/SideWidget.h \
    Core/Exceptions.h \
    Cyclic/CyclicCurves3.h \
    Dependencies/Include/GL/glew.h \
    GA/InterpolatingCurveBiology.h \
    Parametric/ParametricCurves3.h \
    Parametric/ParametricSurfaces3.h \
    Test/TestFunctions.h \
    Test/TorusSurface.h \
    Test/TranguloidTrefoilSurface.h \
    Trigonometric/BicubicBezierPatches.h \
    Trigonometric/TrigonometricBernsteinSurfaces.h \
    Core/GenericCurves3.h \
    Core/HCoordinates3.h \
    Core/Lights.h \
    Core/LinearCombination3.h \
    Core/Materials.h \
    Core/Matrices.h \
    Core/RealSquareMatrices.h \
    Core/TCoordinates4.h \
    Core/TensorProductSurfaces3.h \
    Core/TriangularFaces.h \
    Core/TriangulatedMeshes3.h \
    Core/Constants.h \
    Core/DCoordinates3.h \
    Core/Colors4.h \
    Trigonometric/TrigonometricCurves3.h

SOURCES += \
    GUI/GLWidget.cpp \
    GUI/MainWindow.cpp \
    GUI/SideWidget.cpp \
    main.cpp \
    Cyclic/CyclicCurves3.cpp \
    GA/InterpolatingCurveBiology.cpp \
    Parametric/ParametricCurves3.cpp \
    Parametric/ParametricSurfaces3.cpp \
    Test/TestFunctions.cpp \
    Test/TorusSurface.cpp \
    Test/TranguloidTrefoilSurface.cpp \
    Trigonometric/BicubicBezierPatches.cpp \
    Trigonometric/TrigonometricBernsteinSurfaces.cpp \
    Core/GenericCurves3.cpp \
    Core/Lights.cpp \
    Core/LinearCombination3.cpp \
    Core/Materials.cpp \
    Core/RealSquareMatrices.cpp \
    Core/TensorProductSurfaces3.cpp \
    Core/TriangulatedMeshes3.cpp \
    Trigonometric/TrigonometricCurves3.cpp

CONFIG += console \
          c++11
