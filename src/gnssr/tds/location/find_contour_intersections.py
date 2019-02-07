import numpy as np
from shapely import geometry
import matplotlib.pyplot as plt

def find_contour_intersection(contour1,contour2):
  p1 = contour1.collections[0].get_paths()[0]
  v1 = p1.vertices

  p2 = contour2.collections[0].get_paths()[0]
  v2 = p2.vertices

  poly1 = geometry.LineString(v1)
  poly2 = geometry.LineString(v2)

  intersection = poly1.intersection(poly2)

  return intersection
