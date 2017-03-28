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
	public int id;
	public double x;
	public double y;

	Node(int nid, double nx, double ny) {
		this.id = nid;
		this.x = nx;
		this.y = ny;
	}

	int getID() {
		return id;
	}

	double getX() {
		return x;
	}

	double getY() {
		return y;
	}

	void setID(int nid) {
		this.id = nid;
	}

	void setX(int nx) {
		this.x = nx;
	}

	void setY(int ny) {
		this.y = ny;
	}
}