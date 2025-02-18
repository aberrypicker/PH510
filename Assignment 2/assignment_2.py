# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:49:21 2025.

@author: aaron
"""
import math

# Task 1a

class Vector:
    """
    An object which represents a Vector using coordinates, and can be added,
    subtracted, dot and cross producted, its magnitude determined, and returned
    via print command. Also a conversion scheme for switchign between 
    cartesian and spherical polar coordinates.
    """
    def __init__(self, x, y, z):
        """
        Define the vector by its cartesian coordinates in the x, y,
        and z directions.
        """
        self.x, self.y, self.z = x,y,z
    def __add__(self, other):
        """
        Define the addition of two vectors by adding each component part.
        """
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
    def __sub__(self, other):
        """
        Define the subtraction of two vectors by subtracting their
        component parts.
        """
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)
    def dot(self, other):
        """
        Define the calculation of the dot product between two vectors using
        the dot product formula.
        """
        return self.x * other.x + self.y * other.y + self.z * other.z
    def cross(self, other):
        """
        Define and return the cross product between two vectors.
        """
        return Vector((self.y * other.z - self.z * other.y),
                (self.z * other.x - self.x * other.z),
                (self.x * other.y - self.y *other.x))
    def norm(self):
        """
        Define and return the magnitude of the vector.
        """
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)
    def __repr__(self):
        return f'Vector[x={self.x}, y={self.y}, z={self.z}]'
    # def __eq__(self, other):
    #     print("__eq__ called")
    #     if isinstance(other, tuple):
    #         other = Vector(*other)
    #     if isinstance(other, Vector):
    #         return self.x == other.x and self.y == other.y and self.z == other.z
    #     return NotImplemented
    def conversion(self):
        """
        Converting cartesian coordinates to spherical using the conversion
        equations, plus failsafe for x coord of 0.
        """
        r = math.sqrt(self.x**2 + self.y**2 + self.z**2)

        if self.x == 0:

            theta = math.pi/2
        else:

            theta = math.acos(self.z/r)

        phi = math.atan(self.y/self.x)
        return SphericalPolar(r, theta, phi)


A = Vector(3, 4, 5)
# Task 1b
print()
print("Task 1a & b")
print("Vector A =", A)
# Task 1c
print()
print("Task 1c")
print("The Magnitude of Vector A is", Vector.norm(A))

B = Vector(4, 5, 6)
print("B = ", B)

# Task 1d
print()
print("Task 1d")
C = A + B
print("Vector A + Vector B = [", C, "]")

D = A-B
print("Vector A - Vector B =", D)

# Task 1e
print()
print("Task 1e")
E = Vector.dot(A, B)
print("The Dot Product of A & B is [", E, "]")

F = Vector.cross(A, B)
print("The Cross product of A & B is [", F, "]")


# Task 2
class SphericalPolar(Vector):
    """
    An Object to Represent a Vector in Spherical Coordinates.
    """

    def __init__(self, r, theta, phi):
        """
        Define the Spherical Polar parameters of position, polar angle,
        and azimuthal angle.
        """
        Vector.__init__(self, (r*math.sin(theta)*math.cos(phi)),
                        r*math.sin(theta)*math.sin(phi),
                        r*math.cos(theta))
        self.r, self.theta, self.phi = r, theta, phi
        self.x = self.r*math.sin(self.theta)*math.cos(self.phi)
        self.y = self.r*math.sin(self.theta)*math.sin(self.phi)
        self.z = self.r*math.cos(self.theta)
    def __repr__(self):
        """
        Ensures Spherical Polar coordinate variables are represented when
        printing Polar vectors.
        """
        return f'Spherical Vector[r={self.r:.2f}, theta={self.theta:.2f}, phi={self.phi:.2f}]'


# Task 2a, b: Initialising and Printing Vector in Spherical Polar coords.
P = SphericalPolar(1, 25 * math.pi/180, 43 * math.pi/180)
print()
print("Task 2a + b")
print("Vector P in Spherical Polar Coordinates is", P)

# Task 2c: Printing Magnitude of Spherical Polar Vector.
print()
print("Task 2c")
print("The Magnitude of Vector P is", Vector.norm(P))

# Task 2d: Addition & Subtraction of Spherical Polar Vectors.
Q = SphericalPolar(2, 30 * math.pi/180, 90 * math.pi/180)

R = P + Q
print()
print("Task 2d")
print("Spherical Vector P + Spherical Vector Q =", SphericalPolar.conversion(R))

S = P - Q
print("Spherical Vector P - Spherical Vector Q =", SphericalPolar.conversion(S))

# Task 2e

T = Vector.dot(P, Q)
print()
print("Task 2e")
print("The Dot product between Spherical Polar Vectors P & Q is", T)

U = Vector.cross(P, Q)
print()
print("The Cross product between Spherical Polar Vectors P & Q is", SphericalPolar.conversion(U))
# Task 3a: Triangle Area (Cartesian)
def triangle_area(vertice_1, vertice_2, vertice_3):
    """
    Define a function for returning the area of a Triangle, utilising the cross
    product between two vectors which form joining vertices of the triangle, to
    calculate the area.
    """
    side_1 = vertice_1 - vertice_2
    side_2 = vertice_3 - vertice_2
    return 0.5 * Vector.norm(Vector.cross(side_1, side_2))


A1 = Vector(0, 0, 0)
A2 = Vector(1, 0, 0)
A3 = Vector(0, 1, 0)

Area_A = triangle_area(A1, A2, A3)
print()
print("Task 3a")
print("The Area of Cartesian Triangle A is", Area_A)

B1 = Vector(-1, -1, -1)
B2 = Vector(0, -1, -1)
B3 = Vector(0, 0, -1)

Area_B = triangle_area(B1, B2, B3)
print("The Area of Cartesian Triangle B is", Area_B)

C1 = Vector(1, 0, 0)
C2 = Vector(0, 0, -1)
C3 = Vector(0, 0, 0)

Area_C = triangle_area(C1, C2, C3)
print("The Area of Cartesian Triangle C is", Area_C)

D1 = Vector(0, 0, 0)
D2 = Vector(1, -1, 0)
D3 = Vector(0, 0, 1)

Area_D = triangle_area(D1, D2, D3)
print("The Area of Cartesian Triangle D is", f"{Area_D:.3}")

# Task 3b

def triangle_angle(vertice_1, vertice_2, vertice_3):
    """
    Define a function to use the triangle vertices to return the internal
    angles of the given triangle. The forward and backwards vectors for all
    sides of the triangle are calculated, and then used with dot product
    and magnitude to determine each angle between vectors. The angles, which
    python calculates in radians, are then converted into degrees for 
    ease of visibility.
    """
    side_12 = vertice_2 - vertice_1
    side_23 = vertice_3 - vertice_2
    side_31 = vertice_1 - vertice_3
    side_21 = vertice_1 - vertice_2
    side_32 = vertice_2 - vertice_3
    side_13 = vertice_3 - vertice_1
    angle_1 = math.acos(Vector.dot(side_12, side_13)/
                        (Vector.norm(side_12) * Vector.norm(side_13))) * (180/
                                                                       math.pi)
    angle_2 = math.acos(Vector.dot(side_21, side_23)/
                        (Vector.norm(side_21) * Vector.norm(side_23))) * (180/
                                                                       math.pi)
    angle_3 = math.acos(Vector.dot(side_32, side_31)/
                        (Vector.norm(side_32) * Vector.norm(side_31))) * (180/
                                                                       math.pi)
    return angle_1, angle_2, angle_3


Angles_A = triangle_angle(A1, A2, A3)
print()
print("Task 3b")
print("The Internal Angles of Triangle A are", f"{Angles_A[0]:.1f}",
      f"{Angles_A[1]:.1f}", f"{Angles_A[2]:.1f}")

Angles_B = triangle_angle(B1, B2, B3)
print("The Internal Angles of Triangle B are", f"{Angles_B[0]:.1f}",
      f"{Angles_B[1]:.1f}", f"{Angles_B[2]:.1f}")

Angles_C = triangle_angle(C1, C2, C3)
print("The Internal Angles of Triangle C are", f"{Angles_C[0]:.1f}",
      f"{Angles_C[1]:.1f}", f"{Angles_C[2]:.1f}")

Angles_D = triangle_angle(D1, D2, D3)
print("The Internal Angles of Triangle D are", f"{Angles_D[0]:.1f}",
      f"{Angles_D[1]:.1f}", f"{Angles_D[2]:.1f}")

# Task 3c

P1 = SphericalPolar(0, 0, 0)
P2 = SphericalPolar(1, 0, 0)
P3 = SphericalPolar(1, 90 * math.pi/180, 0)

Area_P = triangle_area(P1, P2, P3)
print()
print("Task 3c")
print("The Area of Spherical Polar Triangle P is", Area_P)
Angles_P = triangle_angle(P1, P2, P3)
print("The Internal Angles of Triangle P are", f"{Angles_P[0]:.1f}",
      f"{Angles_P[1]:.1f}", f"{Angles_P[2]:.1f}")

Q1 = SphericalPolar(1, 0, 0)
Q2 = SphericalPolar(1, 90 * math.pi/180, 0)
Q3 = SphericalPolar(1, 90 * math.pi/180, 180 * math.pi/180)

Area_Q = triangle_area(Q1, Q2, Q3)
print()
print("The Area of Spherical Polar Triangle Q is", f"{Area_Q:.3}")
Angles_Q = triangle_angle(Q1, Q2, Q3)
print("The Internal Angles of Triangle Q are", f"{Angles_Q[0]:.1f}",
      f"{Angles_Q[1]:.1f}", f"{Angles_Q[2]:.1f}")

R1 = SphericalPolar(0, 0, 0)
R2 = SphericalPolar(2, 0, 0)
R3 = SphericalPolar(2, 90 * math.pi/180, 0)

Area_R = triangle_area(R1, R2, R3)
print()
print("The Area of Spherical Polar Triangle R is", f"{Area_R:.3}")
Angles_R = triangle_angle(R1, R2, R3)
print("The Internal Angles of Triangle R are", f"{Angles_R[0]:.1f}",
      f"{Angles_R[1]:.1f}", f"{Angles_R[2]:.1f}")

S1 = SphericalPolar(1, 90 * math.pi/180, 0)
S2 = SphericalPolar(1, 90 * math.pi/180, 180 * math.pi/180)
S3 = SphericalPolar(1, 90 * math.pi/180, 270 * math.pi/180)

Area_S = triangle_area(S1, S2, S3)
print()
print("The Area of Spherical Polar Triangle S is", f"{Area_S:.3}")
Angles_S = triangle_angle(S1, S2, S3)
print("The Internal Angles of Triangle S are", f"{Angles_S[0]:.1f}",
      f"{Angles_S[1]:.1f}", f"{Angles_S[2]:.1f}")
