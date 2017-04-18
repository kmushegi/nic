/*
ACO for TSP - Project 3
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Ant Colony Optimization for the Traveling Salesman Problem,
Project 3. This file is responsible for visualizing the nodes and the tour constructed
by our ACO implementation using the Java Swing framework.

Visualizer currently only works with positive coordinates (+,+)
*/


import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;

public class Visualizer {

	private static ArrayList<Node> nodes;
	private static ArrayList<Node> tour;
	private static int[] bounds;
	private static int scale;
	private static JPanel p;

	private static final int width = 800;
	private static final int height = 600;

    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("ACO TSP Visualizer");
        frame.setLayout(new BorderLayout());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
		initializeNodesAndTour(frame);

        //Display the window.
        frame.setSize(width,height);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    private static void initializeNodesAndTour(JFrame f) {
		p = new JPanel() {
			@Override
			public void paintComponent(Graphics g) {
				super.paintComponent(g);
				Graphics2D g2 = (Graphics2D) g;
				for(int i = 0; i < nodes.size(); i++) {
					NodeDisplay ndt = new NodeDisplay(nodes.get(i).id,nodes.get(i).x,nodes.get(i).y);
					ndt.drawEnv(g2);
				}
				for(int i = 0; i < tour.size(); i++) {
					Node c1 = tour.get(i);
					Node c2 = (i == tour.size()-1) ? tour.get(0) : tour.get(i+1);

					TourLine tlt = new TourLine(c1.x,c1.y,c2.x,c2.y);
					tlt.drawEnv(g2);
				}
			}
		};

		f.getContentPane().add(p, BorderLayout.CENTER);
	}

	//function return format: [minX, maxX, minY, maxY]
	private static int[] findMinAndMaxCoordinates(ArrayList<Node> n) {
		double minX, maxX, minY, maxY;
		minX = minY = Double.MAX_VALUE;
		maxX = maxY = -Double.MAX_VALUE;

		for(int i = 0; i < n.size(); i++) {
			if(n.get(i).x < minX) {
				minX = n.get(i).x;
			}
			if(n.get(i).x > maxX) {
				maxX = n.get(i).x;
			}
			if(n.get(i).y < minY) {
				minY = n.get(i).y;
			}
			if(n.get(i).y > maxY) {
				maxY = n.get(i).y;
			}
		}
		int v[] = new int[4];
		v[0] = (int) Math.round(minX);
		v[1] = (int) Math.round(maxX);
		v[2] = (int) Math.round(minY);
		v[3] = (int) Math.round(maxY);
		return v;
	}

	//re-scale coordinate to fit on canvas.
	private static int scaleCoordinate(double c, double min, double max, double size) {
		double scaledCoordinate;
		scaledCoordinate = c/(max-min) * size;
		return (int)scaledCoordinate;
	}

	private static void convertCoordinates() {
		for(int i = 0; i < nodes.size(); i++) {
			nodes.get(i).x = scaleCoordinate(nodes.get(i).x,bounds[0],bounds[1],width);
			nodes.get(i).y = scaleCoordinate(nodes.get(i).y,bounds[2],bounds[3],height);

			tour.get(i).x = scaleCoordinate(tour.get(i).x,bounds[0],bounds[1],width);
			tour.get(i).y = scaleCoordinate(tour.get(i).y,bounds[2],bounds[3],height);
		}
		tour.get(tour.size()-1).x = scaleCoordinate(tour.get(tour.size()-1).x,bounds[0],bounds[1],width);
		tour.get(tour.size()-1).y = scaleCoordinate(tour.get(tour.size()-1).y,bounds[2],bounds[3],height);
	}

    public static void display(ArrayList<Node> n, ArrayList<Node> t) {
    	nodes = new ArrayList<>(n);
		tour = new ArrayList<>(t);
		bounds = findMinAndMaxCoordinates(nodes);
		scale = Math.max((bounds[1]-bounds[0]),(bounds[3]-bounds[2])); //[minX, maxX, minY, maxY]
		convertCoordinates();

        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGUI();
            }
        });
    }

    static class NodeDisplay {
    	int nid;
    	int xC;
    	int yC;
    	Color nodeColor = Color.black;

    	int w;
    	int h;

    	public NodeDisplay(int id, double x, double y) {
    		nid = id;
    		xC = (int)x;
    		yC = (int)y;
    	}

        public Dimension getPreferredSize() {
            return new Dimension(5, 5);
        }

        protected void drawEnv(Graphics2D g2) {
            g2.setColor(nodeColor);
            g2.setFont(new Font("TimesRoman", Font.PLAIN, 8)); 
            g2.drawString(Integer.toString(nid),xC - 8,yC); 
            g2.fill(new Ellipse2D.Double(xC, yC, 5, 5));
        }
    }

    static class TourLine {
    	int x1;
    	int y1;
    	int x2;
    	int y2;
    	Color nodeColor = Color.red;

    	public TourLine(double nx1, double ny1, double nx2, double ny2) {
    		x1 = (int)nx1;
			y1 = (int)ny1;
			x2 = (int)nx2;
			y2 = (int)ny2;
    	}

    	protected void drawEnv(Graphics2D g2) {
            g2.setColor(nodeColor);
            g2.drawLine(x1,y1,x2,y2);
        }
    }
}