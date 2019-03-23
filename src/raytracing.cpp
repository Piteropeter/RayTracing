#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include "Utilities.h"

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

int height = 1000;
int width = 1000;

float viewport_size = 15.0;
LightSource light_source{3.0, 2.5, 5.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0};
Sphere sphere{1.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.6, 0.7, 0.8, 1.0, 1.0, 1.0, 30.0};
Color global;
float starting_point[3];
float starting_directions[] = {0.0, 0.0, -1.0};
float inter[3];
int inters;
float inters_c[3];
unsigned char pixel[1][1][3];

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
        // if(buffer == "dimensions") {
        //    file >> height >> width;
        //}
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
        //if(step > 0)
        //    ;
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

void Display2(Image& image) {
    int x, y;
    float x_fl, y_fl;
    int im_size_2_x;
    int im_size_2_y;
    im_size_2_x = height / 2;
    im_size_2_y = width / 2;

    Point starting_point2;
    Vector starting_directions2{0, 0, -1};
    Color pixel_to_paint;

    for(y = im_size_2_y; y > -im_size_2_y; y--) {
        for(x = -im_size_2_x; x < im_size_2_x; x++) {
            // std::cout << "X: " << x << "\t\tY: " << y << std::endl;
            x_fl = static_cast<float>(x) / (height / viewport_size);
            y_fl = static_cast<float>(y) / (width / viewport_size);
            // std::cout << "X: " << x_fl << "\t\tY: " << y_fl << std::endl;

            int skonwertowany_x, skonwertowany_y;
            starting_point2.x = x_fl;
            starting_point2.y = y_fl;
            starting_point2.z = viewport_size;
            pixel_to_paint = Trace2(starting_point2, starting_directions2);

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

            // BO BMP JEST ZJEBANE I ZAPISUJE NA ODWRÓT
            image[skonwertowany_x][skonwertowany_y][0] = pixel[0][0][2];
            image[skonwertowany_x][skonwertowany_y][1] = pixel[0][0][1];
            image[skonwertowany_x][skonwertowany_y][2] = pixel[0][0][0];
        }
    }
}

std::array<unsigned char, 14> createBitmapFileHeader(int height, int width, int paddingSize) {
    int fileSize = fileHeaderSize + infoHeaderSize + (bytesPerPixel * width + paddingSize) * height;

    std::array<unsigned char, 14> file_header = {
        0, 0,       /// signature
        0, 0, 0, 0, /// image file size in bytes
        0, 0, 0, 0, /// reserved
        0, 0, 0, 0, /// start of pixel array
    };

    file_header[0] = (unsigned char)('B');
    file_header[1] = (unsigned char)('M');
    file_header[2] = (unsigned char)(fileSize);
    file_header[3] = (unsigned char)(fileSize >> 8);
    file_header[4] = (unsigned char)(fileSize >> 16);
    file_header[5] = (unsigned char)(fileSize >> 24);
    file_header[10] = (unsigned char)(fileHeaderSize + infoHeaderSize);

    return file_header;
}

std::array<unsigned char, 40> createBitmapInfoHeader(int height, int width) {
	std::array<unsigned char, 40> info_header = {
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

    info_header[0] = (unsigned char)(infoHeaderSize);
    info_header[4] = (unsigned char)(width);
    info_header[5] = (unsigned char)(width >> 8);
    info_header[6] = (unsigned char)(width >> 16);
    info_header[7] = (unsigned char)(width >> 24);
    info_header[8] = (unsigned char)(height);
    info_header[9] = (unsigned char)(height >> 8);
    info_header[10] = (unsigned char)(height >> 16);
    info_header[11] = (unsigned char)(height >> 24);
    info_header[12] = (unsigned char)(1);
    info_header[14] = (unsigned char)(bytesPerPixel * 8);

    return info_header;
}

void generateBitmapImage(Image& image, int height, int width, std::string imageFileName) {

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (width * bytesPerPixel) % 4) % 4;

    std::array<unsigned char, 14> fileHeader = createBitmapFileHeader(height, width, paddingSize);
    std::array<unsigned char, 40> infoHeader = createBitmapInfoHeader(height, width);

    FILE* imageFile = fopen(imageFileName.data(), "wb");

    fwrite(fileHeader.data(), 1, fileHeaderSize, imageFile);
    fwrite(infoHeader.data(), 1, infoHeaderSize, imageFile);

    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            fwrite(image[i][j].data(), bytesPerPixel, 1, imageFile);
        }
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);

    // std::ofstream output_file(imageFileName);
    // output_file <<
    // KIEDYŚ
}

void print_help() {
    std::cout << "\nRayTracing Application help:\n\n"
              << "\n\t\t--in-progress--\n\n";
}

std::string extract_name(const std::string& path) {
    std::size_t begin = path.rfind('\\');
    std::size_t end = path.find(".txt");

    return path.substr(begin + 1, end - begin - 1);
}

Parameters get_input_parameters(int args, char* params[]) {
    Parameters parameters;
    // clang-format off
    boost::program_options::options_description desc("Options");
    desc.add_options()
        ("help,h", "Print help messages")
        ("x_resolution,x", boost::program_options::value<std::size_t>(&parameters.x_input_file_resolution), "X resolution of output file, default 1920p")
        ("y_resolution,y", boost::program_options::value<std::size_t>(&parameters.y_input_file_resolution), "Y resolution of output file, default 1080p")
        ("reflections,r", boost::program_options::value<std::size_t>(&parameters.ray_tracing_reflections), "Amount of reflections to calculete for each ray")
        ("input,i", boost::program_options::value<std::string>(&parameters.input_file_path), "Input file")
        ("output,o", boost::program_options::value<std::string>(&parameters.output_file_name), "Output file, default the same name as input");
    // clang-format on
    boost::program_options::positional_options_description pos_desc;
    pos_desc.add("input", 1);

    boost::program_options::variables_map vm;

    boost::program_options::command_line_parser parser{args, params};
    parser.options(desc).positional(pos_desc).allow_unregistered();
    boost::program_options::parsed_options parsed_options = parser.run();
    store(parsed_options, vm);
    notify(vm);

    if(vm.count("help")) {
        parameters.help_mode = true;
    }
    if(vm.count("input")) {
        parameters.is_input_set = true;
    }
    if(vm.count("input") && !vm.count("output")) {
        parameters.output_file_name = extract_name(parameters.input_file_path) + ".bmp";
	}
    return parameters;
}

int main(int args, char* params[]) {
    try {
        auto parameters = get_input_parameters(args, params);
        if(parameters.help_mode || !parameters.is_input_set) {
            print_help();
            return 0;
        } else {
            readFile(parameters.input_file_path);
            Image image;
            image.resize(parameters.x_input_file_resolution);
            for(auto& x : image)
                x.resize(parameters.y_input_file_resolution);

			Display2(image);
            generateBitmapImage(image, height, width, parameters.output_file_name);
            printf("Image generated!!");
            return 0;
        }
    } catch(std::exception& exception) {
        std::cerr << "RayTracing app: " << exception.what() << std::endl;
        return 1;
    }
}