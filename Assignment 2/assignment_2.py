# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:49:21 2025.

@author: aaron
"""
import math

class Vector:
    """
    An object which represents a Vector using coordinates, and can be added,
    subtracted, dot and cross producted, its magnitude determined, and returned
    via print
    """
    
    def __init__(self, x=0, y=0, z=0): #initialises the Vector and its attributes
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
A = Vector(3,4,5)
print("Vector A =",A)
print("The Magnitude of Vector A is", Vector.norm(A))

B = Vector(4,5,6)
print("B = ", B)

C = A + B
print("Vector A + Vector B = [", C, "]")

D = A-B
print("Vector A - Vector B =" , D)

E = Vector.dot(A, B)
print("The Dot Product of A & B is [", E,"]")

F = Vector.cross(A,B)
print("The Cross product of A & B is [", F, "]")



class Spherical_Polar(Vector):
    """An Object to Represent a Vector in Spherical Coordinates"""
    def __init__(self, r, theta, phi):
        """New Parameters"""
        self.r, self.theta, self.phi = r, theta, phi
        self.x = self.r*math.sin(self.theta)*math.cos(self.phi)
        self.y = self.r*math.sin(self.theta)*math.sin(self.phi)
        self.z = self.r*math.cos(self.theta)
    def r(self):
        """Positional Parameter to Spherical Polar Vector"""
        return self.norm()
    def theta(self):
        """Returns the Azimuthal Angle"""
        return math.atan2(self.y, self.x)
    def phi(self):
        """Returns the Polar Angle"""
        return math.acos(self.z/self.norm())
        # if self.degrees == True:
        #     self.theta = self.theta * 180/math.pi
        #     self.phi = self.phi * 180/math.pi

def Triangle_Area(V1, V2, V3):
    Side_1 = V1 - V2
    Side_2 = V3 - V2
    return 0.5 * Vector.norm(Vector.cross(Side_1, Side_2))


#Triangle 1 (Cartesian)

A1 = Vector(0, 0, 0)
A2 = Vector(1, 0, 0)
A3 = Vector(0, 1, 0)

Area_1 = Triangle_Area(A1, A2, A3)
print("The Cartesian Area of Triangle 1 is", Area_1)



