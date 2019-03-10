#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "ArashPatrow/bitmap_image.hpp"
#include "Utilities.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
//#include "Utilities.h"
//#include <vector>

const int bytesPerPixel = 3; /// red, green, blue
const int fileHeaderSize = 14;
const int infoHeaderSize = 40;

typedef float point[3];
using Vector = Point;

Color Phong(const Point& p, const Point& q, const Vector& n, const Vector& d, const int& index, const int step);

double dotProduct(const Point& p1, const Point& p2);

Vector Normalization(const Point& p);
Vector Reflect(const Point& p, const Point& q, const Vector& d, const Vector& n);
Point viewer_v{0.0, 0.0, 1.0}; // wektor kierunku obserwacji
int im_size_x = 5000;
int im_size_y = 5000; // Nadpisywane wczytaniem z pliku

const int height = 5000;
const int width = 5000;
std::array<std::array<std::array<unsigned char, bytesPerPixel>, width>, height> image;

float viewport_size = 15.0;
LightSource light_source{3.0, 2.5, 5.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0};
Sphere sphere{1.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.6, 0.7, 0.8, 1.0, 1.0, 1.0, 30.0};
Color global;
float starting_point[3];                        // punkt, z którego wychodzi promieñ
float starting_directions[] = {0.0, 0.0, -1.0}; // wektor opisuj¹cy kierunek promienia
float inter[3];                                 // wspó³rzêdne (x,y,z) punktu przeciêcia promienia i sfery
int inters;                                     // zmienna okreœlaj¹ca, czy sfera zosta³a przeciêta przez
float inters_c[3];                              // sk³adowe koloru dla oœwietlonego punktu na powierzchni sfery
unsigned char pixel[1][1][3];                   // sk³adowe koloru rysowanego piksela

const int MAX_STEPS = 3;

std::vector<Sphere> spheres;
std::vector<LightSource> lights;
Color background;

void readFile(const std::string& fileName) {
    std::ifstream file;
    file.open(fileName);
    std::string buffer = "";
    while(!file.eof()) {
        file >> buffer;
        if(buffer == "dimensions") {
            file >> im_size_x >> im_size_y;
        }
        if(buffer == "background") {
            file >> background.r >> background.g >> background.b;
        }
        if(buffer == "global") {
            file >> global.r >> global.g >> global.b;
        }
        if(buffer == "sphere") {
            Sphere newSphere{};
            file >> newSphere.r;
            file >> newSphere.x0 >> newSphere.y0 >> newSphere.z0;
            file >> newSphere.KsR >> newSphere.KsG >> newSphere.KsB; // specular
            file >> newSphere.KdR >> newSphere.KdG >> newSphere.KdB; // diffuse
            file >> newSphere.KaR >> newSphere.KaG >> newSphere.KaB; // ambient
            file >> newSphere.n;

            spheres.push_back(newSphere);
        }
        if(buffer == "source") {
            LightSource newLight{};
            file >> newLight.xs >> newLight.ys >> newLight.zs;
            file >> newLight.IsR >> newLight.IsG >> newLight.IsB;
            file >> newLight.IdR >> newLight.IdG >> newLight.IdB;
            file >> newLight.IaR >> newLight.IaG >> newLight.IaB;

            lights.push_back(newLight);
        }
    }
    file.close();
}

Point Intersect(const Point& p, const Vector& v, std::pair<int, unsigned int>& status) {
    Point intersection{0, 0, 0};
    double a, b, c, d, r, lenght = DBL_MAX;
    for(auto i = 0u; i < lights.size(); i++) {
        a = lights[i].xs - p.x;
        b = lights[i].ys - p.y;
        c = lights[i].zs - p.z;

        if((a / v.x) == (b / v.y) && (b / v.y) == (c / v.z)) {
            intersection.x = lights[i].xs;
            intersection.y = lights[i].ys;
            intersection.z = lights[i].zs;
            status.first = 2;
            status.second = i;
        }
    }
    if(!status.first) {
        for(auto i = 0u; i < spheres.size(); i++) {
            a = v.x * v.x + v.y * v.y + v.z * v.z;
            b = 2 * ((p.x - spheres[i].x0) * v.x + (p.y - spheres[i].y0) * v.y + (p.z - spheres[i].z0) * v.z);
            c = (p.x - spheres[i].x0) * (p.x - spheres[i].x0) + (p.y - spheres[i].y0) * (p.y - spheres[i].y0) +
                (p.z - spheres[i].z0) * (p.z - spheres[i].z0) - spheres[i].r * spheres[i].r;

            d = b * b - 4 * a * c;
            if(d >= 0) // jest co najmniej jeden punkt przeciêcia
            {
                r = (-b - sqrt(d)) / (2 * a); // parametr dla bli¿szego punktu przeciêcia
                if(r > 0 && r < lenght) {
                    intersection.x = p.x + r * v.x; // wspó³rzêdne punktu przeciêcia
                    intersection.y = p.y + r * v.y;
                    intersection.z = p.z + r * v.z;
                    lenght = sqrt(pow(intersection.x - p.x, 2.0) + pow(intersection.y - p.y, 2.0) + pow(intersection.z - p.z, 2.0));

                    status.first = 1;
                    status.second = i;
                    // jest punkt przeciêcia
                }
            }
        }
    }
    if(status.first) {
        // std::cout << "\nX: " << intersection.x << "\tY: " << intersection.y << "\tZ: " << intersection.z;
        return intersection;
    }
    // promieñ nie przecina sfery
    status.first = 0;
    status.second = 0;
    return intersection;
}

Vector Normal(const Point& q, unsigned int i) {
    double x = 0, y = 0, z = 0;
    x = (q.x - spheres[i].x0) / spheres[i].r;
    y = (q.y - spheres[i].y0) / spheres[i].r;
    z = (q.z - spheres[i].z0) / spheres[i].r;
    return Vector{x, y, z};
}

Color Trace2(Point p, Vector d, int step = 0) {
    std::pair<int, unsigned int> data{0, 0}; // first - rodzaj outputu (0 - nic, 1 - sfera, 2 - œwiat³o); second - index
    // auto& status = data.first;
    if(step > MAX_STEPS) // przeanalizowano ju¿ zadan¹ liczbê poziomów drzewa
        // return background;
        return Color{0, 0, 0};

    Color local, reflected;          // sk³adowe koloru
    Point q;                         // wspó³rzêdne punktu
    Vector n, r;                     // wspó³rzêdne wektora
    q = Intersect(p, d, data);       // obliczenie pubnktu przeciêcia promienia i obiektu sceny
    if(data.first == 0 && step == 0) // nic nie zosta³o trafione
    {
        return background;
    }
    if(data.first == 1) {
        if(step > 0)
            ;
        n = Normal(q, data.second);                   // obliczenie wektora normalnego w punkcie q
        r = Reflect(p, q, d, n);                      // obliczenie wektoru odbicia promienia w punkcie q
        local = Phong(p, q, n, d, data.second, step); // oblieczenie oœwietlenia lokalnego w punkcie q
        reflected = Trace2(q, r, step + 1);           // obliczenie "reszty" oœwietlenia dla punktu q
        return (local + reflected);                   // obliczenie ca³kowitego oœwietlenia dla q
    }
    return Color{0, 0, 0};
    // return background;
}

Vector Reflect(const Point& p, const Point& q, const Vector& d, const Vector& n) {
    double x = 0, y = 0, z = 0;
    // Vector ref = { -d.x, -d.y, -d.z };
    Vector ref = {p.x - q.x, p.y - q.y, p.z - q.z};
    ref = Normalization(ref);
    double n_dot_1 = dotProduct(ref, n);
    x = 2 * n.x * n_dot_1 - ref.x;
    y = 2 * n.y * n_dot_1 - ref.y;
    z = 2 * n.z * n_dot_1 - ref.z;
    if(dotProduct(Vector{x, y, z}, Vector{x, y, z}) > 0)
        return Normalization(Vector{x, y, z});
    return Vector{x, y, z};
}

Color Phong(const Point& p, const Point& q, const Vector& n, const Vector& d, const int& index, const int step) {
    double x = 0, y = 0, z = 0;
    Vector ref;
    ref = {-d.x, -d.y, -d.z};
    // ref = { 0, 0, 1 };
    ref = Normalization(ref);

    const double a = 1, b = 0.01, c = 0.001;
    // sprawdzenie czy punkt na powierzchni sfery jest oœwietlany przez Ÿród³o
    for(auto i = 0u; i < lights.size(); i++) {
        Vector light_vector{lights[i].xs - q.x, lights[i].ys - q.y, lights[i].zs - q.z};
        const double d2 = sqrt(dotProduct(light_vector, light_vector));
        light_vector = Normalization(light_vector);

        double n_dot_l = dotProduct(light_vector, n);
        Vector reflection_vector{2 * n_dot_l * n.x - light_vector.x, 2 * n_dot_l * n.y - light_vector.y, 2 * n_dot_l * n.z - light_vector.z};
        reflection_vector = Normalization(reflection_vector);

        double v_dot_r = dotProduct(reflection_vector, ref);
        if(v_dot_r < 0) // obserwator nie widzi oœwietlanego punktu
            v_dot_r = 0;

        if(n_dot_l > 0) // punkt jest oœwietlany, oœwietlenie wyliczane jest ze wzorów dla modelu Phonga
        {
            x = x +
                (spheres[index].KdR * lights[i].IdR * n_dot_l + spheres[index].KsR * lights[i].IsR * pow(double(v_dot_r), spheres[index].n)) *
                    (1.00 / (a + b * d2 + c * d2 * d2)) +
                spheres[index].KaR * lights[i].IaR;
            y = y +
                (spheres[index].KdG * lights[i].IdG * n_dot_l + spheres[index].KsG * lights[i].IsG * pow(double(v_dot_r), spheres[index].n)) *
                    (1.00 / (a + b * d2 + c * d2 * d2)) +
                spheres[index].KaG * lights[i].IaG;
            z = z +
                (spheres[index].KdB * lights[i].IdB * n_dot_l + spheres[index].KsB * lights[i].IsB * pow(double(v_dot_r), spheres[index].n)) *
                    (1.00 / (a + b * d2 + c * d2 * d2)) +
                spheres[index].KaB * lights[i].IaB;
        }
        // punkt nie jest oœwietlany, uwzglêdniane jest tylko œwiat³o rozproszone
        x = x + spheres[index].KaR * global.r;
        y = y + spheres[index].KaG * global.g;
        z = z + spheres[index].KaB * global.b;
    }
    return Color{x, y, z};
}

Vector Normalization(const Point& p) {
    double x = 0, y = 0, z = 0;
    auto d = 0.0;
    // d = p.x * p.x + p.y * p.y + p.z * p.z;
    d = dotProduct(p, p);
    d = sqrt(d);
    if(d > 0.0) {
        x = p.x / d;
        y = p.y / d;
        z = p.z / d;
        return Vector{x, y, z};
    } else {
        return Vector{p.x, p.y, p.z};
    }
}

double dotProduct(const Point& p1, const Point& p2) {
    return (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
}

void Display2(void) {
    int x, y;         // pozycja rysowanego piksela "ca³kowitoliczbowa"
    float x_fl, y_fl; // pozycja rysowanego piksela "zmiennoprzecinkowa"
    int im_size_2_x;  // po³owa rozmiaru obrazu w pikselach
    int im_size_2_y;
    im_size_2_x = im_size_x / 2; // obliczenie po³owy rozmiaru obrazu w pikselach
    im_size_2_y = im_size_y / 2;
    // glClear(GL_COLOR_BUFFER_BIT);
    // glFlush();
    Point starting_point2;
    Vector starting_directions2{0, 0, -1};
    Color pixel_to_paint;

    // rysowanie pikseli od lewego górnego naro¿nika do prawego dolnego naro¿nika
    for(y = im_size_2_y; y > -im_size_2_y; y--) {
        for(x = -im_size_2_x; x < im_size_2_x; x++) {
            //std::cout << "X: " << x << "\t\tY: " << y << std::endl;
            x_fl = static_cast<float>(x) / (im_size_x / viewport_size);
            y_fl = static_cast<float>(y) / (im_size_y / viewport_size);
            //std::cout << "X: " << x_fl << "\t\tY: " << y_fl << std::endl;

            int skonwertowany_x, skonwertowany_y;
            // przeliczenie pozycji(x,y) w pikselach na pozycjê "zmiennoprzecinkow¹" w oknie obserwatora
            starting_point2.x = x_fl;
            starting_point2.y = y_fl;
            starting_point2.z = viewport_size;
            // wyznaczenie pocz¹tku œledzonego promienia dla rysowanego piksela
            pixel_to_paint = Trace2(starting_point2, starting_directions2);
            // obliczenie punktu przeciêcia ze sfer¹
            // glRasterPos3f(x_fl, y_fl, 0);
            // inkrementacja pozycji rastrowej dla rysowania piksela

			skonwertowany_x = x + im_size_2_x;
			skonwertowany_y = y + im_size_2_y - 1;


            if(pixel_to_paint.r > 1)
                pixel[0][0][0] = 255;
            else
                pixel[0][0][0] = pixel_to_paint.r * 255;

            if(pixel_to_paint.g > 1)
                pixel[0][0][1] = 255;
            else
                pixel[0][0][1] = pixel_to_paint.g * 255;

            if(pixel_to_paint.b > 1)
                pixel[0][0][2] = 255;
            else
                pixel[0][0][2] = pixel_to_paint.b * 255;
			

			image[skonwertowany_x][skonwertowany_y][0] = pixel[0][0][0];
			image[skonwertowany_x][skonwertowany_y][1] = pixel[0][0][1];
			image[skonwertowany_x][skonwertowany_y][2] = pixel[0][0][2];
            // glDrawPixels(1, 1, GL_RGB, GL_UNSIGNED_BYTE, pixel);
            // narysowanie kolejnego piksela na ekranie
            // glFlush();
        }
    }
}

// void Myinit(void) {
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    glOrtho(-viewport_size / 2, viewport_size / 2, -viewport_size / 2, viewport_size / 2, -viewport_size / 2, viewport_size / 2);
//    glMatrixMode(GL_MODELVIEW);
//}

// int main(void) {
//    readFile("scene.txt");
//    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
//    glutInitWindowSize(im_size_x, im_size_y);
//    glutCreateWindow("Ray Casting 2000");
//    Myinit();
//    glutDisplayFunc(Display2);
//    glutMainLoop();
//    return 0;
//}

// void print_help() {
//    std::cout << "\nPEA Project Application help:\n\n"
//              << "Synopsis:\n"
//              << "\tapplication [OPTIONS]... [INPUT_FILE]\n\n"
//              << "Description:\n"
//              << "\t-h, --help\n\t\thelp screen\n"
//              << "\t-d, --diversification\n\t\tenables diversifiaction of starting path for tabu search\n"
//              << "\t-s [SEED], --seed [SEED]\n\t\tsets seed for random number generator"
//              << "\n\t\tdefault seed is 2018\n"
//              << "\t-l [SIZE], --tabu_list [SIZE]\n\t\tsize of tabu list"
//              << "\n\t\tdefault size is 100\n"
//              << "\t-p [SIZE], --population [SIZE]\n\t\tsize of population in genetic algorithm"
//              << "\n\t\tdefault size is 1000\n"
//              << "\t-t [SECONDS], --time [SECONDS]\n\t\tduration of one tabu search cycle"
//              << "\n\t\tdefault time is one second\n"
//              << "\t-b [ITERATIONS], --benchmark [ITERATIONS]\n\t\tinstead results application returns time measurments"
//              << "\n\t\titerations specify how many measurments will be taken\n"
//              << "\t-i [ITERATIONS], --iterations [ITERATIONS]\n\t\titerations specify how many iterations will be made in one cycle in tabu search"
//              << "\n\t\tdefault amount is 100\n"
//              << "\t-c [PROBABILITY], --crossing [PROBABILITY]\n\t\tprobability of crossing between two paths in genetic algorithm"
//              << "\n\t\tdefault chance is 0.9\n"
//              << "\t-n [PROBABILITY], --mutation [PROBABILITY]\n\t\tprobability of path mutation in genetic algorithm"
//              << "\n\t\tdefault chance is 0.01\n"
//              << "\t-o [MODE], --output [MODE]\n\t\toutput mode"
//              << "\n\t\tavailable modes are console (default mode) and file (file is generated)\n"
//              << "\t-m [MODE], --mode [MODE]\n\t\tsolving mode"
//              << "\n\t\tavailable modes are brute (brute force), bnb (branch and bound), dynamic (dynamic programming), tabu (tabu search), genetic (genetic algorithm)"
//              << "\n\t\tdefault mode is tabu search\n\n";
//}

// void print_results(const std::string& results, const Reader::Parameters& params, const std::string& file_name) {
//    if(params.file_mode) {
//        std::string mode;
//        switch(params.solving_mode) {
//        case Reader::MODE::BRUTE_FORCE:
//            mode = "_brute";
//            break;
//        case Reader::MODE::BRANCH_AND_BOUND:
//            mode = "_bnb";
//            break;
//        case Reader::MODE::DYNAMIC:
//            mode = "_dynamic";
//            break;
//        case Reader::MODE::TABU:
//            mode = "_tabu" + params.benchmark_repetitions;
//            break;
//        case Reader::MODE::GENETIC:
//            mode = "_genetic" + params.benchmark_repetitions;
//        default:
//            throw std::runtime_error("Solving mode not set!");
//        }
//        std::string output;
//        if(params.benchmark_mode)
//            output = file_name + mode + "_" + std::to_string(params.benchmark_repetitions) + "_output.txt";
//        else
//            output = file_name + mode + "_output.txt";
//
//        std::ofstream output_file(output);
//        output_file << results;
//        output_file.close();
//    } else
//        std::cout << results;
//}
//
// void proceed(const Reader::Parameters& parameters) {
//    Reader::Graph graph = Reader::read_graph(parameters.path);
//
//    auto results = Solver::solve(graph, parameters);
//
//    print_results(results, parameters, graph.name);
//}
//
// Reader::Parameters process_input(int args, char* params[]) {
//    Reader::Parameters output_parameters;
//
//    std::string solving_mode = "";
//    std::string output_mode = "";
//
//    boost::program_options::options_description desc("Options");
//    desc.add_options()
//        ("help,h", "Print help messages")
//        ("mode,m", boost::program_options::value<std::string>(&solving_mode), "Choose solving mode")
//        ("benchmark,b", boost::program_options::value<std::size_t>(&output_parameters.benchmark_repetitions), "Benchmark mode")
//        ("output,o", boost::program_options::value<std::string>(&output_mode), "Output mode")
//        ("iterations,i", boost::program_options::value<std::size_t>(&output_parameters.algorithm_iterations), "Iterations for tabu search or genetic algorithm")
//        ("diversification,d", "Enables diversification for tabu search")
//        ("seed,s", boost::program_options::value<std::size_t>(&output_parameters.seed), "Seed for random number generator")
//        ("tabu_list,l", boost::program_options::value<std::size_t>(&output_parameters.tabu_list_size), "Tabu list size")
//		("population,p", boost::program_options::value<std::size_t>(&output_parameters.population_size), "Population size for genetic algorithm")
//        ("time,t", boost::program_options::value<double>(&output_parameters.time), "Time constraint for tabu search or genetic algorithm")
//        ("crossing,c", boost::program_options::value<double>(&output_parameters.crossing_chance), "Crossing chance for genetic algorithm")
//        ("mutation,n", boost::program_options::value<double>(&output_parameters.mutation_chance), "Mutation chance for genetic algorithm")
//        ("file", boost::program_options::value<std::string>(&output_parameters.path), "Input file");
//
//    boost::program_options::positional_options_description pos_desc;
//    pos_desc.add("file", 1);
//
//    boost::program_options::variables_map vm;
//
//    boost::program_options::command_line_parser parser{args, params};
//    parser.options(desc).positional(pos_desc).allow_unregistered();
//    boost::program_options::parsed_options parsed_options = parser.run();
//    store(parsed_options, vm);
//    notify(vm);
//
//    if(vm.count("help")) {
//        output_parameters.help_mode = true;
//    }
//    if(vm.count("diversification")) {
//        output_parameters.diversification = true;
//    }
//    if(vm.count("file")) {
//        output_parameters.is_path_set = true;
//    }
//    if(vm.count("benchmark")) {
//        output_parameters.benchmark_mode = true;
//    }
//    if(vm.count("mode")) {
//        if(solving_mode == "brute")
//            output_parameters.solving_mode = Reader::MODE::BRUTE_FORCE;
//        else if(solving_mode == "bnb")
//            output_parameters.solving_mode = Reader::MODE::BRANCH_AND_BOUND;
//        else if(solving_mode == "dynamic")
//            output_parameters.solving_mode = Reader::MODE::DYNAMIC;
//        else if(solving_mode == "tabu")
//            output_parameters.solving_mode = Reader::MODE::TABU;
//        else if(solving_mode == "genetic")
//            output_parameters.solving_mode = Reader::MODE::GENETIC;
//        else
//            throw std::logic_error("invalid solving mode!");
//    }
//    if(vm.count("output")) {
//        if(output_mode == "file")
//            output_parameters.file_mode = true;
//        else if(output_mode == "console")
//            output_parameters.file_mode = false;
//        else
//            throw std::logic_error("invalid output mode!");
//    }
//    return output_parameters;
//}

// int main(int args, char* params[]) {
//    /*try {
//        Reader::Parameters parameters = process_input(args, params);
//
//        if(parameters.help_mode) {
//            print_help();
//            return 0;
//        }
//        if(parameters.is_path_set)
//            proceed(parameters);
//        else
//            std::cout << "input file is not set! \nIf you need help run application with --help flag." << std::endl;
//
//    } catch(std::exception& exception) {
//        std::cerr << "Application: " << exception.what() << std::endl;
//        return 1;
//    }*/
//
//
//	int w = 128;
//    int h = 128;
//
//
//
//
//	//int* red = new int[w][h];
//    //int* green = new int[w][h];
//    //int* blue = new int[w][h];
//    std::array<std::array<std::uint8_t, 128>, 128> red;
//    std::array<std::array<std::uint8_t, 128>, 128> green;
//    std::array<std::array<std::uint8_t, 128>, 128> blue;
//
//	for(int i = 0; i < w; i++)
//        for(int j = 0; j < h; j++) {
//
//            red[i][j] = 0xaa;
//            green[i][j] = 0xaa;
//            blue[i][j] = 0xaa;
//				}
//
//
//	FILE* f;
//    unsigned char* img = NULL;
//    int filesize = 54 + 3 * w * h; // w is your image width, h is image height, both int
//
//    img = (unsigned char*)malloc(3 * w * h);
//    memset(img, 0, 3 * w * h);
//
//    for(int i = 0; i < w; i++) {
//        for(int j = 0; j < h; j++) {
//            int x = i;
//            int y = (h - 1) - j;
//            int r = red[i][j] * 255;
//            int g = green[i][j] * 255;
//            int b = blue[i][j] * 255;
//            if(r > 255)
//                r = 255;
//            if(g > 255)
//                g = 255;
//            if(b > 255)
//                b = 255;
//            img[(x + y * w) * 3 + 2] = (unsigned char)(r);
//            img[(x + y * w) * 3 + 1] = (unsigned char)(g);
//            img[(x + y * w) * 3 + 0] = (unsigned char)(b);
//        }
//    }
//
//    unsigned char bmpfileheader[14] = {'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0};
//    unsigned char bmpinfoheader[40] = {40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0};
//    unsigned char bmppad[3] = {0, 0, 0};
//
//    bmpfileheader[2] = (unsigned char)(filesize);
//    bmpfileheader[3] = (unsigned char)(filesize >> 8);
//    bmpfileheader[4] = (unsigned char)(filesize >> 16);
//    bmpfileheader[5] = (unsigned char)(filesize >> 24);
//
//    bmpinfoheader[4] = (unsigned char)(w);
//    bmpinfoheader[5] = (unsigned char)(w >> 8);
//    bmpinfoheader[6] = (unsigned char)(w >> 16);
//    bmpinfoheader[7] = (unsigned char)(w >> 24);
//    bmpinfoheader[8] = (unsigned char)(h);
//    bmpinfoheader[9] = (unsigned char)(h >> 8);
//    bmpinfoheader[10] = (unsigned char)(h >> 16);
//    bmpinfoheader[11] = (unsigned char)(h >> 24);
//
//    f = fopen("img.bmp", "wb");
//    fwrite(bmpfileheader, 1, 14, f);
//    fwrite(bmpinfoheader, 1, 40, f);
//    for(int i = 0; i < h; i++) {
//        fwrite(img + (w * (h - i - 1) * 3), 3, w, f);
//        fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
//    }
//
//    free(img);
//    fclose(f);
//
//
//
//
//
//
//
//
//
//    return 0;
//}
//

void generateBitmapImage(unsigned char* image, int height, int width, char* imageFileName);
unsigned char* createBitmapFileHeader(int height, int width, int paddingSize);
unsigned char* createBitmapInfoHeader(int height, int width);

int main() {
    readFile("scene.txt");
    // unsigned char image[height][width][bytesPerPixel];
    char* imageFileName = "bitmapImage.bmp";

    //int i, j;
    //for(i = 0; i < height; i++) {
    //    for(j = 0; j < width; j++) {
    //        image[i][j][2] = (unsigned char)((double)i / height * 255);                 /// red
    //        image[i][j][1] = (unsigned char)((double)j / width * 255);                  /// green
    //        image[i][j][0] = (unsigned char)(((double)i + j) / (height + width) * 255); /// blue
    //    }
    //}
    Display2();
    generateBitmapImage((unsigned char*)image.data(), height, width, imageFileName);
    printf("Image generated!!");
}

void generateBitmapImage(unsigned char* image, int height, int width, char* imageFileName) {

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (width * bytesPerPixel) % 4) % 4;

    unsigned char* fileHeader = createBitmapFileHeader(height, width, paddingSize);
    unsigned char* infoHeader = createBitmapInfoHeader(height, width);

    FILE* imageFile = fopen(imageFileName, "wb");

    fwrite(fileHeader, 1, fileHeaderSize, imageFile);
    fwrite(infoHeader, 1, infoHeaderSize, imageFile);

    int i;
    for(i = 0; i < height; i++) {
        fwrite(image + (i * width * bytesPerPixel), bytesPerPixel, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);
    // free(fileHeader);
    // free(infoHeader);
}

unsigned char* createBitmapFileHeader(int height, int width, int paddingSize) {
    int fileSize = fileHeaderSize + infoHeaderSize + (bytesPerPixel * width + paddingSize) * height;

    static unsigned char fileHeader[] = {
        0, 0,       /// signature
        0, 0, 0, 0, /// image file size in bytes
        0, 0, 0, 0, /// reserved
        0, 0, 0, 0, /// start of pixel array
    };

    fileHeader[0] = (unsigned char)('B');
    fileHeader[1] = (unsigned char)('M');
    fileHeader[2] = (unsigned char)(fileSize);
    fileHeader[3] = (unsigned char)(fileSize >> 8);
    fileHeader[4] = (unsigned char)(fileSize >> 16);
    fileHeader[5] = (unsigned char)(fileSize >> 24);
    fileHeader[10] = (unsigned char)(fileHeaderSize + infoHeaderSize);

    return fileHeader;
}

unsigned char* createBitmapInfoHeader(int height, int width) {
    static unsigned char infoHeader[] = {
        0, 0, 0, 0, /// header size
        0, 0, 0, 0, /// image width
        0, 0, 0, 0, /// image height
        0, 0,       /// number of color planes
        0, 0,       /// bits per pixel
        0, 0, 0, 0, /// compression
        0, 0, 0, 0, /// image size
        0, 0, 0, 0, /// horizontal resolution
        0, 0, 0, 0, /// vertical resolution
        0, 0, 0, 0, /// colors in color table
        0, 0, 0, 0, /// important color count
    };

    infoHeader[0] = (unsigned char)(infoHeaderSize);
    infoHeader[4] = (unsigned char)(width);
    infoHeader[5] = (unsigned char)(width >> 8);
    infoHeader[6] = (unsigned char)(width >> 16);
    infoHeader[7] = (unsigned char)(width >> 24);
    infoHeader[8] = (unsigned char)(height);
    infoHeader[9] = (unsigned char)(height >> 8);
    infoHeader[10] = (unsigned char)(height >> 16);
    infoHeader[11] = (unsigned char)(height >> 24);
    infoHeader[12] = (unsigned char)(1);
    infoHeader[14] = (unsigned char)(bytesPerPixel * 8);

    return infoHeader;
}