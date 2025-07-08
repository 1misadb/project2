#include <GL/glew.h>
#include <GL/glut.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <fstream>

#include "OpenglCallbacks.h"
#include "PreProcessing.h"
#include "algorithms.h"

// Global state replicated from original GUI version
GLdouble vertices[64][6];
int vertexIndex = 0;
int tx, ty, tw, th;
int handle;
vector<Polygon_2> (*nfps)(vector<Polygon_2>, vector<Polygon_2>);
char problemName[1024] = "";
Layout layout;
int saveAll = 0;
int staticHeuristic = WIDER;
list<Point_2> piecesOrdered;
map<Point_2, GLuint> layoutNFPsTest;
map<Point_2, map<Point_2, DrawingWithRotations>> drawingNFPsTest;
map<Point_2, DrawingWithRotations> drawingPolysTest;
int iteration = 0;
vector<int> piecesAvaliability;
int stockList = 0;
GLfloat *blue = nullptr;
Polygon_with_holes_2 **polygons = nullptr;
int numberOfPolygons = 0;
int nrOfRotations = 0;
DrawingWithRotations currentDrawingPolys;
map<Point_2, DrawingWithRotations> currentDrawingNFPs;
std::ofstream testfile;

static void initGL() {
    glShadeModel(GL_SMOOTH);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_CULL_FACE);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glClearColor(0, 0, 0, 0);
    glClearStencil(0);
    glClearDepth(1.0f);
    glDepthFunc(GL_LEQUAL);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    std::cout << "==== CLI Nesting started ====\n";
    std::cout << "Input file: " << argv[1] << "\n";
    std::cout << "Output file: " << argv[2] << "\n";

    std::strncpy(problemName, argv[1], sizeof(problemName) - 1);
    problemName[sizeof(problemName) - 1] = '\0';
    const char* outputFile = argv[2];

    // default configuration
    saveAll = 1;               // use dynamic algorithm
    staticHeuristic = WIDER;   // default static heuristic
    nfps = sortEdgesWithDecomposition;

    std::cout << "[1] Initializing shared memory...\n";
    initSharedMem();

    std::cout << "[2] Initializing OpenGL context...\n";
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);
    glutInitWindowSize(500, 200);
    glutCreateWindow(argv[0]);

    if (GLEW_OK != glewInit()) {
        std::cerr << "Failed to init GLEW" << std::endl;
        return 1;
    }
    initGL();

    std::cout << "[3] Loading layout...\n";
    std::string inFile(problemName);
    std::transform(inFile.begin(), inFile.end(), inFile.begin(), ::tolower);
    if (inFile.rfind(".dat") != std::string::npos)
        layout = loadConfigurationFileIrregularProblem(problemName);
    else
        layout = loadLayout(problemName);

    std::cout << "[4] Preprocessing...\n";
    preProcessing();

    std::cout << "[5] Building layout...\n";
    if (saveAll) {
        std::cout << "[5.1] Using dynamic algorithm\n";
        buildDynamicLayout(outputFile);
    } else {
        std::cout << "[5.2] Using static algorithm\n";
        buildStaticLayout(outputFile);
    }

    std::cout << "[6] Clearing shared memory...\n";
    clearSharedMem();

    std::cout << "==== Nesting complete. Output saved to: " << outputFile << " ====\n";
    return 0;
}
