import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.io.*;

public class Visualizer {

	private static ArrayList<Node> nodes;//= new ArrayList<>();
	private static ArrayList<Node> tour;//= new ArrayList<>();
	JPanel p;

    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("ACO TSP Visualizer");
        frame.setLayout(new BorderLayout());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
		initializeNodes(frame);

        //Display the window.
        frame.setSize(1000,800);
        // frame.pack();
        frame.setVisible(true);
    }

    private static void initializeNodes(JFrame f) {
		JPanel p = new JPanel() {
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

    private static void processProblemLine(String[] tokens) {
		try {
			Integer.parseInt(tokens[0]); //make sure its a node line
			tokens = Utility.formatNodeInput(tokens);

			Node temp = new Node(Integer.parseInt(tokens[0]),
							Double.parseDouble(tokens[1]),
							Double.parseDouble(tokens[2]));
			nodes.add(temp);
		} catch (NumberFormatException e) {
			// System.err.format("NumberFormatException %s%n",e);
		}
	}

    private static void readProblem(String fp) {
		try(BufferedReader br = new BufferedReader(new FileReader(fp))) {
			String line;
			while((line = br.readLine()) != null) {
				String[] tokens = line.trim().split(" "); //trim ws & tokenize
				processProblemLine(tokens);
			}
		} catch (IOException e) {
			System.err.format("IOException %s%n",e);
			Logger.printErrorAndExit("Problem File Not Found. Exiting");
		}
	}

    public static void display(ArrayList<Node> n, ArrayList<Node> t) {
    	// readProblem(args[0]);
    	nodes = new ArrayList<>(n);
		tour = new ArrayList<>(t);

		// for(int i = 0; i < tour.size(); i++) {
		// 	System.out.println(tour.get(i).x + " " + tour.get(i).y);
		// }

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

    	public NodeDisplay(int id, double x, double y) {
    		nid = id;
    		xC = (int)x;
    		yC = (int)y;
    	}

        public Dimension getPreferredSize() {
            return new Dimension(5, 5);
        }

        protected void drawEnv(Graphics2D g2) {
            // super.paintComponent(g2);
            g2.setColor(nodeColor);
            g2.setFont(new Font("TimesRoman", Font.PLAIN, 8)); 
            g2.drawString(Integer.toString(nid),xC/2 - 5,yC/2-2); 
            g2.fill(new Ellipse2D.Double(xC/2, yC/2, 5, 5));
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
            // super.paintComponent(g2);
            g2.setColor(nodeColor);
            g2.drawLine(x1/2,y1/2,x2/2,y2/2);
        }
    }
}