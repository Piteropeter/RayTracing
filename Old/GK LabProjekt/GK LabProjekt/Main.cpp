/*************************************************************************************/
// Program rysuje obraz barwnej sfery oœwietlonej przez jedno barwne Ÿród³o œwiat³a.
// Œrodek sfery znajduje siê w œrodku uk³adu wspó³rzêdnych.
// Do obliczania oœwietlenia wykorzystywany jest model Phonga.
/*************************************************************************************/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <string>
#include <iostream>
#include <fstream>
#include "Objects.h"
#include <vector>
typedef float point[3];
using Vector = Point;
/*************************************************************************************/
// Prototypy u¿ywanych funkcji
/*************************************************************************************/
// Funkcja obliczaj¹ca punkt przeciêcia promienia i powierzchni sfery
/*************************************************************************************/
// Funkcja obliczaj¹ca oœwietlenie punktu na powierzchni sfery wed³ug modelu Phonga
Color Phong(const Point& p, const Point& q, const Vector& n, const Vector& d, const int& index, const int step);
/*************************************************************************************/
// Funkcja obliczaj¹ca iloczyn skalarny dwóch wektorów
double dotProduct(const Point& p1, const Point& p2);
/*************************************************************************************/
// Funkcja normalizuj¹ca wektor
Vector Normalization(const Point& p);
/*************************************************************************************/
Vector Reflect(const Point& p, const Point& q, const Vector& d, const Vector& n);
// Zmienne globalne
Point viewer_v{ 0.0, 0.0, 1.0 }; // wektor kierunku obserwacji
/*************************************************************************************/
// Rozmiar obrazu w pikselach (obraz jest kwadratem)
int im_size_x = 1000;
int im_size_y = 1000; //Nadpisywane wczytaniem z pliku
/*************************************************************************************/
// Rozmiar okna obserwatora
float viewport_size = 15.0;
//float viewport_size = 3.0;
/*************************************************************************************/
// Po³o¿enie i parametry Ÿród³a œwiat³a
LightSource light_source{ 3.0, 2.5, 5.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
/*************************************************************************************/
// Promieñ i parametry rysowanej sfery
//float sphere_radius = 1.0;
//float sphere_specular[] = {0.8f, 0.8f, 0.8f};
//float sphere_diffuse[] = {0.6f, 0.7f, 0.8f};
//float sphere_ambient[] = {1.0, 1.0, 1.0};
//float sphere_specularhininess = 30.0;
Sphere sphere{ 1.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.6, 0.7, 0.8, 1.0, 1.0, 1.0, 30.0 };
/*************************************************************************************/
// Parametry œwiat³a rozproszonego
Color global;
/*************************************************************************************/
// Parametry "œledzonego" promienia
float starting_point[3]; // punkt, z którego wychodzi promieñ
float starting_directions[] = {0.0, 0.0, -1.0}; // wektor opisuj¹cy kierunek promienia
/*************************************************************************************/
// Zmienne pomocnicze
float inter[3]; // wspó³rzêdne (x,y,z) punktu przeciêcia promienia i sfery
int inters; // zmienna okreœlaj¹ca, czy sfera zosta³a przeciêta przez
float inters_c[3]; // sk³adowe koloru dla oœwietlonego punktu na powierzchni sfery
GLubyte pixel[1][1][3]; // sk³adowe koloru rysowanego piksela

const int MAX_STEPS = 3;

std::vector<Sphere> spheres;
std::vector<LightSource> lights;
Color background;

void readFile(const std::string& fileName) {
    std::ifstream file;
    file.open(fileName);
    std::string buffer = "";
    while (!file.eof()) {
        file >> buffer;
        if (buffer == "dimensions") {
           file >> im_size_x >> im_size_y;
        }
        if(buffer == "background") {
            file >> background.r >> background.g >> background.b;
        }
        if (buffer == "global") {
            file >> global.r >> global.g >> global.b;
        }
        if (buffer == "sphere") {
			Sphere newSphere{};
			file >> newSphere.r;
            file >> newSphere.x0 >> newSphere.y0 >> newSphere.z0;
            file >> newSphere.KsR >> newSphere.KsG >> newSphere.KsB; //specular
            file >> newSphere.KdR >> newSphere.KdG >> newSphere.KdB; //diffuse
            file >> newSphere.KaR >> newSphere.KaG >> newSphere.KaB; //ambient
            file >> newSphere.n;

			spheres.push_back(newSphere);
        }
        if (buffer == "source") {
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

/*************************************************************************************/
// Funkcja oblicza punkt przeciêcia promienia i powierzchni sfery
// Argument p jest punktem pocz¹tkowym promienia a v wektorem opisuj¹cym
// kierunek biegu promienia
// Funkcja zwraca 1 jeœli promieñ przecina sferê, 0 gdy nie przecina.
/*************************************************************************************/
Point Intersect(const Point& p, const Vector& v, std::pair<int, unsigned int>& status) {
	Point intersection{0,0,0};
	double a, b, c, d, r, lenght = DBL_MAX;
	for (auto i = 0u; i < lights.size(); i++) {
		a = lights[i].xs - p.x;
		b = lights[i].ys - p.y;
		c = lights[i].zs - p.z;

		if ((a / v.x) == (b / v.y) && (b / v.y) == (c / v.z)) {
			intersection.x = lights[i].xs;
			intersection.y = lights[i].ys;
			intersection.z = lights[i].zs;
			status.first = 2;
			status.second = i;
		}
	}
	if (!status.first) {
		for (auto i = 0u; i < spheres.size(); i++) {
			a = v.x * v.x + v.y * v.y + v.z * v.z;
			b = 2 * ((p.x - spheres[i].x0) * v.x + (p.y - spheres[i].y0) * v.y + (p.z - spheres[i].z0) * v.z);
			c = (p.x - spheres[i].x0) * (p.x - spheres[i].x0) + (p.y - spheres[i].y0) * (p.y - spheres[i].y0) + (p.z - spheres[i].z0) * (p.z - spheres[i].z0) - spheres[i].r * spheres[i].r;
			/*c = dotProduct(p, p) - 2 * (spheres[i].x0 * p.x + spheres[i].y0 * p.y + spheres[i].z0 * p.z) +
				dotProduct(Point{ spheres[i].x0, spheres[i].y0, spheres[i].z0 }, Point{ spheres[i].x0, spheres[i].y0, spheres[i].z0 }) -
				pow(spheres[i].r, 2.0);*/

			d = b * b - 4 * a * c;
			if (d >= 0) // jest co najmniej jeden punkt przeciêcia
			{
				r = (-b - sqrt(d)) / (2 * a); // parametr dla bli¿szego punktu przeciêcia
				if (r > 0 && r < lenght) {
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
	if (status.first) {
		//std::cout << "\nX: " << intersection.x << "\tY: " << intersection.y << "\tZ: " << intersection.z;
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
	return Vector{ x, y, z };
}

Color Trace2(Point p, Vector d, int step = 0)
{
	std::pair<int, unsigned int> data{0, 0};				  //first - rodzaj outputu (0 - nic, 1 - sfera, 2 - œwiat³o); second - index
	//auto& status = data.first;
	if (step > MAX_STEPS)                     //przeanalizowano ju¿ zadan¹ liczbê poziomów drzewa
		//return background;
		return Color{ 0, 0, 0 };
    
    Color local, reflected;                    //sk³adowe koloru
    Point q;                                    //wspó³rzêdne punktu
	Vector n, r;                                //wspó³rzêdne wektora
    q = Intersect(p, d, data);          //obliczenie pubnktu przeciêcia promienia i obiektu sceny

	//if (data.first == 2) //trafione Ÿród³o œwiat³a
		//return Color{lights[data.second].IsR, lights[data.second].IsG, lights[data.second].IsB}; //never happened
		//throw std::exception("IT WORKS");
	if (data.first == 0 && step == 0) //nic nie zosta³o trafione
	{
		return background;
	}
	if (data.first == 1) {
		if (step > 0);
		n = Normal(q, data.second);                           //obliczenie wektora normalnego w punkcie q
		r = Reflect(p, q, d, n);                     //obliczenie wektoru odbicia promienia w punkcie q
		local = Phong(p, q, n, d, data.second, step);  //oblieczenie oœwietlenia lokalnego w punkcie q
		reflected = Trace2(q, r, step + 1);      //obliczenie "reszty" oœwietlenia dla punktu q
		return (local + reflected);              //obliczenie ca³kowitego oœwietlenia dla q
	}
	return Color{ 0, 0, 0 };
	//return background;

} 

Vector Reflect(const Point& p, const Point& q, const Vector& d, const Vector& n)
{
	double x = 0, y = 0, z = 0;
	//Vector ref = { -d.x, -d.y, -d.z };
	Vector ref = { p.x - q.x, p.y - q.y, p.z - q.z };
	ref = Normalization(ref);
	double n_dot_1 = dotProduct(ref, n);
	x = 2 * n.x * n_dot_1 - ref.x;
	y = 2 * n.y * n_dot_1 - ref.y;
	z = 2 * n.z * n_dot_1 - ref.z;
	if(dotProduct(Vector{ x, y, z }, Vector{ x, y, z }) > 0)
		return Normalization(Vector{x, y, z});
	return Vector{ x, y, z };

}

/*************************************************************************************/
// Funkcja oblicza oœwietlenie punktu na powierzchni sfery u¿ywaj¹c modelu Phonga
/*************************************************************************************/
Color Phong(const Point& p, const Point& q, const Vector& n , const Vector& d, const int& index, const int step) {
	double x = 0, y = 0, z = 0;
	Vector ref;
	ref = { -d.x, -d.y, -d.z };
	//ref = { 0, 0, 1 };
	ref = Normalization(ref);

	const double a = 1, b = 0.01, c = 0.001;
	// sprawdzenie czy punkt na powierzchni sfery jest oœwietlany przez Ÿród³o
	for (auto i = 0u; i < lights.size(); i++) {
		Vector light_vector{lights[i].xs - q.x, lights[i].ys - q.y, lights[i].zs - q.z };
		const double d2 = sqrt(dotProduct(light_vector, light_vector));
		light_vector = Normalization(light_vector);

		double n_dot_l = dotProduct(light_vector, n);
		Vector reflection_vector{
			2 * n_dot_l * n.x - light_vector.x, 
			2 * n_dot_l * n.y - light_vector.y,
			2 * n_dot_l * n.z - light_vector.z
		};
		reflection_vector = Normalization(reflection_vector);

		double v_dot_r = dotProduct(reflection_vector, ref);
		if (v_dot_r < 0) // obserwator nie widzi oœwietlanego punktu
			v_dot_r = 0;

		if (n_dot_l > 0) // punkt jest oœwietlany, oœwietlenie wyliczane jest ze wzorów dla modelu Phonga
		{
			x = x + (spheres[index].KdR * lights[i].IdR * n_dot_l + spheres[index].KsR * lights[i].IsR * pow(double(v_dot_r), spheres[index].n)) * (1.00 / (a + b * d2 + c * d2 * d2))
				+ spheres[index].KaR * lights[i].IaR;
			y = y + (spheres[index].KdG * lights[i].IdG * n_dot_l + spheres[index].KsG * lights[i].IsG * pow(double(v_dot_r), spheres[index].n)) * (1.00 / (a + b * d2 + c * d2 * d2))
				+ spheres[index].KaG * lights[i].IaG;
			z = z + (spheres[index].KdB * lights[i].IdB * n_dot_l + spheres[index].KsB * lights[i].IsB * pow(double(v_dot_r), spheres[index].n)) * (1.00 / (a + b * d2 + c * d2 * d2))
				+ spheres[index].KaB * lights[i].IaB;
		}
		// punkt nie jest oœwietlany, uwzglêdniane jest tylko œwiat³o rozproszone  
		x = x + spheres[index].KaR * global.r;
		y = y + spheres[index].KaG * global.g;
		z = z + spheres[index].KaB * global.b;		
	}
	return Color{x, y, z};
}

/*************************************************************************************/
// Funkcja przeprowadza normalizacjê wektora
/*************************************************************************************/

Vector Normalization(const Point& p) {
	double x = 0, y = 0, z = 0;
	auto d = 0.0;
	//d = p.x * p.x + p.y * p.y + p.z * p.z;
	d = dotProduct(p, p);
	d = sqrt(d);
	if (d > 0.0) {
		x = p.x / d;
		y = p.y / d;
		z = p.z / d;
		return Vector{ x, y, z };
	}
	else {
		return Vector{ p.x, p.y, p.z };
	}
}

/*************************************************************************************/
// Funkcja oblicza iloczyn skalarny wektorów
/*************************************************************************************/
double dotProduct(const Point& p1, const Point& p2) {
    return (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
}

/*************************************************************************************/
// Funkcja rysuj¹ca obraz oœwietlonej sceny
/*************************************************************************************/
void Display2(void) {
	int x, y; // pozycja rysowanego piksela "ca³kowitoliczbowa"
	float x_fl, y_fl; // pozycja rysowanego piksela "zmiennoprzecinkowa"
	int im_size_2_x; // po³owa rozmiaru obrazu w pikselach
	int im_size_2_y;
	im_size_2_x = im_size_x / 2; // obliczenie po³owy rozmiaru obrazu w pikselach
	im_size_2_y = im_size_y / 2;
	glClear(GL_COLOR_BUFFER_BIT);
	glFlush();
	Point starting_point2;
	Vector starting_directions2{ 0, 0, -1 };
	Color pixel_to_paint;
	// rysowanie pikseli od lewego górnego naro¿nika do prawego dolnego naro¿nika
	for (y = im_size_2_y; y > -im_size_2_y; y--) {
		for (x = -im_size_2_x; x < im_size_2_x; x++) {
			x_fl = static_cast<float>(x) / (im_size_x / viewport_size);
			y_fl = static_cast<float>(y) / (im_size_y / viewport_size);
			// przeliczenie pozycji(x,y) w pikselach na pozycjê "zmiennoprzecinkow¹" w oknie obserwatora
			starting_point2.x = x_fl;
			starting_point2.y = y_fl;
			starting_point2.z = viewport_size;
			// wyznaczenie pocz¹tku œledzonego promienia dla rysowanego piksela
			pixel_to_paint = Trace2(starting_point2, starting_directions2);
			// obliczenie punktu przeciêcia ze sfer¹
			glRasterPos3f(x_fl, y_fl, 0);
			// inkrementacja pozycji rastrowej dla rysowania piksela
			if (pixel_to_paint.r > 1)
				pixel[0][0][0] = 255;
			else 
				pixel[0][0][0] = pixel_to_paint.r * 255;

			if (pixel_to_paint.g > 1)
				pixel[0][0][1] = 255;
			else
				pixel[0][0][1] = pixel_to_paint.g * 255;

			if (pixel_to_paint.b > 1)
				pixel[0][0][2] = 255;
			else
				pixel[0][0][2] = pixel_to_paint.b * 255;

			glDrawPixels(1, 1, GL_RGB, GL_UNSIGNED_BYTE, pixel);
			// narysowanie kolejnego piksela na ekranie
			glFlush();
		}
	}
}
/*************************************************************************************/
// Funkcja inicjalizuj¹ca definiuj¹ca sposób rzutowania
/*************************************************************************************/
void Myinit(void) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-viewport_size / 2, viewport_size / 2, -viewport_size / 2, viewport_size / 2, -viewport_size / 2, viewport_size / 2);
    glMatrixMode(GL_MODELVIEW);
}
/*************************************************************************************/
// Funkcja g³ówna
/*************************************************************************************/
int main(void) {
    readFile("scene.txt");
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(im_size_x, im_size_y);
    glutCreateWindow("Ray Casting 2000");
    Myinit();
    glutDisplayFunc(Display2);
    glutMainLoop();
    return 0;
}
/*************************************************************************************/