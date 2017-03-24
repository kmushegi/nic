/*
ACO for TSP - Project 3
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

The code in this file contains the implementation of the Node class.
*/

import java.io.*;
import java.util.*;

class Node {
	public double x;
	public double y;

	Node(int nx, int ny) {
		this.x = nx;
		this.y = ny;
	}
}