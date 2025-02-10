# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:49:21 2025

@author: aaron
"""
import math

class Vector:
    """An object which represents a Vector using coordinates"""
    def __init__(self, x=0, y=0, z=0): #initialises the Vector and its attributes
        self.x, self.y, self.z = x,y,z
    def __add__(self, other): #addition of two vectors
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
    def __sub__(self, other): #subtraction of two vectors
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)
    def __dot__(self, other): #dot product between two vectors
        return self.x * other.x + self.y * other.y + self.z * other.z
    def __cross__(self, other): #cross product between two vectors
        return Vector((self.y * other.z - self.z * other.y),
                (self.z * other.x - self.x * other.z),
                (self.x * other.y - self.y *other.x))
    def __abs__(self): #magnitude of vector
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

class Sphere(Vector):
    """An Object to Represent a Vector in Spherical Coordinates"""
A = Vector(3,4,5)
print("Vector A = [", A.x, A.y, A.z, "]")
print("The Magnitude of Vector A is", Vector.__abs__(A))

B = Vector(4,5,6)
print("B = [", B.x, B.y, B.z, "]")

C = Vector.__add__(A, B)
ctest = A+B
print("Vector A + Vector B = [", C, "]")

D = Vector.__sub__(A,B)
print("Vector A - Vector B = [D" , D,"]")

E = Vector.__dot__(A, B)
print("The Dot Product of A & B is [", E,"]")

F = Vector.__cross__(A,B)
print("The Cross product of A & B is [", F, "]")
