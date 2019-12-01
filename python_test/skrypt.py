from math import sin, cos

staly_tekst1 = """background 0.3 0.3 0.3
global 0.1   0.1    0.1
sphere 0.7   3.0    0.0  -5.0      0.8 0.2 0.0  0.7 1.0 0.0  0.2 0.1 0.2  40
sphere 0.7  -3.0    0.0  -5.0      0.8 0.2 0.0  0.7 1.0 0.0  0.2 0.1 0.2  40
sphere 2.0   0.0    0.0  -1.0      0.8 0.1 0.0  0.8 0.1 0.0  0.2 0.1 0.2  40
sphere 2.0   0.0   -5.0  -3.0      0.8 0.2 0.0  0.0 0.7 1.0  0.2 0.1 0.2  40
sphere 2.0   0.0    5.0  -3.0      0.8 0.2 0.0  0.0 0.7 1.0  0.2 0.1 0.2  40
sphere 2.0  -5.0    2.5  -3.0      0.8 0.2 0.0  0.0 0.7 1.0  0.2 0.1 0.2  40
sphere 2.0  -5.0   -2.5  -3.0      0.8 0.2 0.0  0.0 0.7 1.0  0.2 0.1 0.2  40
sphere 2.0   5.0   -2.5  -3.0      0.8 0.2 0.0  0.0 0.7 1.0  0.2 0.1 0.2  40
sphere 2.0   5.0    2.5  -3.0      0.8 0.2 0.0  0.0 0.7 1.0  0.2 0.1 0.2  40
source       0.0    0.0  15.0      0.2 0.2 0.2  0.4 0.4 0.4  0.2 0.2 0.2
source      -5.0    0.0  10.0      0.2 0.2 0.2  1.0 0.0 1.0  0.3 0.3 0.1
source       5.0    0.0  10.0      0.2 0.2 0.2  1.0 0.0 1.0  0.3 0.3 0.1
source       5.0    0.0  12.0      0.2 0.2 0.2  0.0 1.0 1.0  0.4 0.5 0.3
source      -5.0    0.0  12.0      0.2 0.2 0.2  0.0 1.0 1.0  0.4 0.5 0.3"""

staly_tekst2 = """background 0.0 0.0 0.0
global 0.01 0.01 0.01
sphere 7.5   0.0    0.0  0.0      0.99 0.99 0.99  0.99 0.99 0.99  0.99 0.99 0.99  1000000"""
#source       0.0    0.0  10.0      1.0 1.0 1.0  1.0 1.0 1.0  1.0 1.0 1.0"""

# source      -5.0    0.0  10.0      0.2 0.2 0.2  1.0 0.0 1.0  0.3 0.3 0.1
# source       5.0    0.0  10.0      0.2 0.2 0.2  1.0 0.0 1.0  0.3 0.3 0.1
# source       5.0    0.0  12.0      0.2 0.2 0.2  0.0 1.0 1.0  0.4 0.5 0.3
# source      -5.0    0.0  12.0      0.2 0.2 0.2  0.0 1.0 1.0  0.4 0.5 0.3

class Sphere:
    def __init__(self, r, x, y, z, lx1, ly1, lz1, lx2, ly2, lz2, lx3, ly3, lz3, n, level, offset):
        self.r = r
        self.x = x
        self.y = y
        self.z = z
        self.lx1 = lx1
        self.ly1 = ly1
        self.lz1 = lz1
        self.lx2 = lx2
        self.ly2 = ly2
        self.lz2 = lz2
        self.lx3 = lx3
        self.ly3 = ly3
        self.lz3 = lz3
        self.n = n
        self.level = level
        self.offset = offset

class Light:
    def __init__(self, x, y, z, lx1, ly1, lz1, lx2, ly2, lz2, lx3, ly3, lz3):
        self.x = x
        self.y = y
        self.z = z
        self.lx1 = lx1
        self.ly1 = ly1
        self.lz1 = lz1
        self.lx2 = lx2
        self.ly2 = ly2
        self.lz2 = lz2
        self.lx3 = lx3
        self.ly3 = ly3
        self.lz3 = lz3


offset1 = [0.10,
           0.15,
           0.09,
           0.16,
           0.08,
           0.07,
           0.20,
           0.06,
           0.11,
           0.12,
           0.13]
 
files_needed = int(81.00 / min(offset1))

print("WITAM W TURBO SKRYPCIE!")

spheres1 = {Sphere(5.5, -20.0, 0.0, -10.0, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 40, 0, offset1[0]),
           Sphere(2.5, -20.0, 3.0, -7.0, 0.9, 0.9, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 50, 0, offset1[1]),
           Sphere(0.5, -20.0, 2.0, -4.0, 0.9, 0.0, 0.9, 0.1, 0.5, 0.1, 0.1, 0.1, 0.1, 50, 0, offset1[2]),
           Sphere(1.5, -20.0, -1.0, 0.0, 0.9, 0.0, 0.5, 0.1, 0.2, 0.1, 0.1, 0.3, 0.1, 50, 0, offset1[3]),
           Sphere(4.5, -20.0, -5.0, -8.0, 0.1, 0.3, 0.9, 0.1, 0.5, 0.9, 0.1, 0.1, 0.1, 50, 0, offset1[4]),
           Sphere(0.5, -20.0, 0.5, 1.0, 0.1, 0.3, 0.9, 0.1, 0.5, 0.9, 0.1, 0.1, 0.1, 50, 0, offset1[5]),
           Sphere(2.5, -20.0, -5.0, -12.0, 0.2, 0.3, 0.9, 0.3, 0.5, 0.0, 0.2, 0.1, 0.1, 80, 0, offset1[6]),
           Sphere(1.5, -20.0, -1.0, -15.0, 0.1, 0.2, 0.3, 0.1, 0.5, 0.7, 0.1, 0.3, 0.4, 50, 0, offset1[7]),
           Sphere(5.5, -20.0, -4.5, -42.0, 0.4, 0.1, 0.9, 0.2, 0.1, 0.6, 0.5, 0.1, 0.4, 60, 0, offset1[8]),
           Sphere(1.5, -20.0, 2.0, -8.0, 0.1, 0.3, 0.4, 0.1, 0.5, 0.9, 0.1, 0.4, 0.1, 50, 0, offset1[9]),
           Sphere(1.5, -20.0, 3.0, -18.0, 0.5, 0.3, 0.4, 0.4, 0.5, 0.5, 0.4, 0.1, 0.1, 40, 0, offset1[10])}

spheres2 = {}

lights = {Light(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0),
        #   Light(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0),
        #   Light(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0),
          Light(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0)}


def write_sphere(f, sphere):
    f.write(
        "\nsphere\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d" % (sphere.r, sphere.x, sphere.y, sphere.z,
                                                                              sphere.lx1, sphere.ly1, sphere.lz1,
                                                                              sphere.lx2, sphere.ly2, sphere.lz2,
                                                                              sphere.lx3, sphere.ly3, sphere.lz3,
                                                                              sphere.n))
      
def write_circulating_light(f, light, i, offset):
    f.write(
        "\nsource\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (light.x, light.y + sin(float(i)/30) * 1000000 * offset, light.z + cos(float(i)/30) * 1000000 * offset,
                                                                        light.lx1, light.ly1, light.lz1,
                                                                        light.lx2, light.ly2, light.lz2,
                                                                        light.lx3, light.ly3, light.lz3))

# for i in range(files_needed):
offset = 1
for i in range(1515):
    f = open("output/%04d.txt" % (i + 1), "w+")
    f.write(staly_tekst2)
    for light in lights:
        write_circulating_light(f, light, i, offset)
        offset = offset * -1

    # for sphere in spheres1:
    #     if sphere.level != 2:
    #         if sphere.x >= 20.0:
    #             sphere.level = 1
    #         if sphere.x <= -21.0:
    #             sphere.level = 2
    #         if sphere.level == 0:
    #             sphere.x = sphere.x + sphere.offset
    #         if sphere.level == 1:
    #             sphere.x = sphere.x - sphere.offset
    #     write_sphere(f, sphere)

print("WYGENEROWANE!")
