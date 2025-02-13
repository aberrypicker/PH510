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
    via print
    """
    def __init__(self, x, y, z): #initialises the Vector and its attributes
        self.x, self.y, self.z = x,y,z
    def __add__(self, other):
        """For addition of two vectors"""
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
    def __sub__(self, other):
        """ For subtraction of two vectors"""
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)
    def dot(self, other):
        """To calculate the dot product between two vectors"""
        return self.x * other.x + self.y * other.y + self.z * other.z
    def cross(self, other):
        """To determine the cross product between two vectors"""
        return Vector((self.y * other.z - self.z * other.y),
                (self.z * other.x - self.x * other.z),
                (self.x * other.y - self.y *other.x))
    def norm(self):
        """Determining and returning the magnitude of the vector"""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)
    def __repr__(self):
        return f'Vector[x={self.x}, y={self.y}, z={self.z}]'
    def __eq__(self, other):
        print("__eq__ called")
        if isinstance(other, tuple):
            other = Vector(*other)
        if isinstance(other, Vector):
            return self.x == other.x and self.y == other.y and self.z == other.z
        return NotImplemented


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
print("Vector A - Vector B =" , D)

# Task 1e
print()
print("Task 1e")
E = Vector.dot(A, B)
print("The Dot Product of A & B is [", E,"]")

F = Vector.cross(A,B)
print("The Cross product of A & B is [", F, "]")


# Task 2
class SphericalPolar(Vector):
    """An Object to Represent a Vector in Spherical Coordinates."""

    def __init__(self, r, theta, phi):
        """Define the Spherical Polar parameters."""
        self.r, self.theta, self.phi = r, theta, phi
        self.x = self.r*math.sin(self.theta)*math.cos(self.phi)
        self.y = self.r*math.sin(self.theta)*math.sin(self.phi)
        self.z = self.r*math.cos(self.theta)
        # if self.degrees == True:
        #     self.theta = self.theta * 180/math.pi
        #     self.phi = self.phi * 180/math.pi

Q = SphericalPolar(1, 45, 30)



# Task 3a: Triangle Area (Cartesian)
def triangle_area(vertice_1, vertice_2, vertice_3):
    """
    Function for determining area of a Triangle
    using paralellogram cross product
    """
    side_1 = vertice_1 - vertice_2
    side_2 = vertice_3 - vertice_2
    return 0.5 * Vector.norm(Vector.cross(side_1, side_2))


A1 = Vector(0, 0, 0)
A2 = Vector(1, 0, 0)
A3 = Vector(0, 1, 0)

Area_1 = triangle_area(A1, A2, A3)
print()
print("Task 3a")
print("The Cartesian Area of Triangle A is", Area_1)

B1 = Vector(-1, -1, -1)
B2 = Vector(0, -1, -1)
B3 = Vector(0, 0, -1)

Area_B = triangle_area(B1, B2, B3)
print("The Cartesian Area of Triangle B is", Area_B)

C1 = Vector(1, 0, 0)
C2 = Vector(0, 0, -1)
C3 = Vector(0, 0, 0)

Area_C = triangle_area(C1, C2, C3)
print("The Cartesian Area of Triangle C is", Area_C)

D1 = Vector(0, 0, 0)
D2 = Vector(1, -1, 0)
D3 = Vector(0, 0, 1)

Area_D = triangle_area(D1, D2, D3)
print("The Cartesian Area of Triangle D is", f"{Area_D:.3f}")

# Task 3b

def triangle_angle(vertice_1, vertice_2, vertice_3):
    """
    Function to use the triangle vertices to determine the internal
    angles of the given triangle.
    """
    side_1 = vertice_1 - vertice_2
    side_2 = vertice_3 - vertice_2
    side_3 = vertice_2 + vertice_3 
    angle_1 = math.acos(Vector.dot(side_1, side_2)/ (Vector.norm(side_1) * Vector.norm(side_2))) * (180/math.pi)
    angle_2 = math.acos(Vector.dot(side_3, side_2)/
                        (Vector.norm(side_3) * Vector.norm(side_2))) * (180/math.pi)
    angle_3 = 180 - angle_1 - angle_2
    return angle_1, angle_2, angle_3

Angles_A = triangle_angle(A1, A2, A3)
print()
print("Task 3b")
print("The Internal Angles of Triangle A are", Angles_A)

Angles_B = triangle_angle(B1, B2, B3)
print("The Internal Angles of Triangle B are", Angles_B)

Angles_C = triangle_angle(C1, C2, C3)
print("The Internal Angles of Triangle C are", Angles_C)

Angles_D = triangle_angle(D1, D2, D3)
print("The Internal Angles of Triangle D are", Angles_D)

# Task 3c