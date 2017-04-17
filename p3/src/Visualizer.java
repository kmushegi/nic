import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.io.*;

public class Visualizer {

	private static ArrayList<Node> nodes = new ArrayList<>();
	JPanel p;

    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("ACO TSP Visualizer");
        frame.setLayout(new BorderLayout());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        // JPanel panel = new JPanel();
        // panel.setPreferredSize(new Dimension(800,600));

		// for(int i = 0; i < 5; i++) {
		// 	NodeDisplay nt = new NodeDisplay((int)nodes.get(i).x,(int)nodes.get(i).y);
		// 	frame.getContentPane().add(nt);
		// }
		initializeNodes(frame);
		// panel.add(nt);

        //Display the window.
        frame.setSize(800,600);
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

    public static void main(String[] args) {
    	readProblem(args[0]);

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
            return new Dimension(15, 15);
        }

        protected void drawEnv(Graphics2D g2) {
            // super.paintComponent(g2);
            g2.setColor(nodeColor);
            g2.fill(new Ellipse2D.Double(xC/2, yC/2, 15, 15));
        }
    }
}