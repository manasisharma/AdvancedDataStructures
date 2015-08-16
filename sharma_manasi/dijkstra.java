/* author: manasi */

import java.util.Scanner;
import java.util.Set;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.HashSet;

/* main class */
public class dijkstra 
{ 
	static int n, nEdges = 0;   /* n is the number of nodes, nEdges is the number of edges */
	
	public static int sourceNode = 0;
	
	
	public static void main(String[] args) throws FileNotFoundException
	{
		Vertex[] g;
		Edge[] edges = null;
		int n, j;
               		double b;
	
		if(args[0].equals("-r"))			//random mode
                {
                   
                    int r1 = Integer.parseInt(args[1]);
                    double r2 = Double.parseDouble(args[2]);
                    int r3 = Integer.parseInt(args[3]);
	
/*  randomly generating graph with r1 nodes and r2 density in order to calulate the shortest path from the source vertex r3 */
                    random(r1,r2,r3);
                }
		
		else
                {
                       
                        String filename = args[1];
/* read the graph information into the destination */
                        edges = readFile(filename);
                    
                    
                    n = getNumberOfNodes(edges);
/* generate the graph with n and edges as number of nodes and edges */
                    g = generateGraph(edges, n);
                    int[][] edgeCost = new int[n][n];
                    
/* calculating time for the shortest path from the source vertex using simple scheme*/
                    if(args[0].equals("-s"))		
                    {
                       
                        b = System.currentTimeMillis();
                        edgeCost = simpleScheme(g);

                        b = System.currentTimeMillis() - b;

 //time calculation for simple scheme
                      
                        
                    }

/* calculating time for the shortest path from the source vertex using fibonacci heap scheme*/
                   
	 else if(args[0].equals("-f"))	
                    {
                      
                        b = System.currentTimeMillis();
                        edgeCost = fheapScheme(g);
                        b = System.currentTimeMillis() - b; 

//time calculation for fibonacci scheme
                    
                    }
                    else{
                   
 /* for arguments other than -s and - f */
                   
	 System.out.println("Invalid Arguments"");
                    System.exit(0);
                   
                    }
              
 /*displaying the overall cost from the source node*/     
                        
	for(j=0; j<n; j++)
                        {
                            if(edgeCost[sourceNode][j] == Integer.MAX_VALUE)
                                System.out.print("?\t");
                            else
                                System.out.println(edgeCost[sourceNode][j]);// + "\t" + "// Cost from "+sourceNode+" to "+j);
                        }
                        System.out.println();
                  
                }
	}
	
/*mode : random - - generate a graph randomly - - and  compare the performance of two different schemes namely : simple and fibonacci heap */

	public static void random(int a, double y, int z)
	{
		n = a; 								

//number of node in graph
 
		double density = (y * 0.01); 

/* y is the density receievd as an argument to the function random() */

		Edge[] edges;
		Edge tempEdge;
		Vertex[] g = null;
		int i, start, end;
		double perc;
		Random gen = new Random();
		boolean isConnected = false,duplicateEdge;
		double t1,t2;
		System.out.println("Number of vertices\t" + "Density\t\t");
		perc = (density * 100.0);
		System.out.println(n + "\t\t\t" + perc +"%");
		n= n/5;
		nEdges = (int) (n * (n-1) * density * 0.5);	
		isConnected = false;    

 /*initializing isConnected variable to false in order to mark the graph disconnected and then initializing it to true once we ensure the graph is connected */
				
							
				edges = new Edge[nEdges];
				while(!isConnected)			

//execute the two schemes only when the graph is a connected

				{
					System.out.println("Generating graph... ");
					i = 0;
					while(i < nEdges)		

//create n*(n-1)*density edges

					{
						start = gen.nextInt(n);
						end = gen.nextInt(n);
						{
							tempEdge = new Edge(start, end, gen.nextInt(1000) + 1);
							duplicateEdge = false;
							for(Edge e : edges)
							{      
								if(tempEdge.equals(e))	

 //checking for duplicates : check whether the new generated edge already exists

								{	
									duplicateEdge = true;
									break;
								}
							}
							
							if(duplicateEdge == false)
							{
								edges[i] = tempEdge;
								i ++;
							}

						}
					}

/*generating the graph with n nodes and edges*/

					g = generateGraph(edges, n);
					isConnected = isConnected(g);		

//test whether the generated graph is connected
				
					}
				
				t1 = System.currentTimeMillis();
				simpleScheme(g);
				
				t1 = System.currentTimeMillis() - t1;
				t2 = System.currentTimeMillis();
				fheapScheme(g);

				t2 = System.currentTimeMillis() - t2;
				
				
/* displaying the final time taken by the two schemes */
                                
				
				System.out.println("Time taken by Simple scheme "+ t1 +" Milliseconds");
				System.out.println("Time taken by Fibonacci heap Scheme "+ t2 + " Milliseconds");
	}
	
	
/* readFile () reads the desired edge information from a file wherein  filename is the file that stores the edge information and returns the edges also throws FileNotFoundException */
	
	public static Edge[] readFile(String filename) throws FileNotFoundException
	{
		Scanner in = new Scanner(new File(filename));

//the above scanner class is getting the string from values
		
		int first;
		first = in.nextInt();
		sourceNode = first;
		int count = 0;
		ArrayList<Integer> U1 = new ArrayList<Integer>(); //node1
		ArrayList<Integer> U2 = new ArrayList<Integer>(); //node2
		ArrayList<Integer> U3 = new ArrayList<Integer>(); //weight

		while (in.hasNextInt()) {
			if (count == 0) { 

//scanner is reading the first line in the file
 
				n = in.nextInt();
				nEdges = in.nextInt();
				count++;
			}
			U1.add(in.nextInt());
			U2.add(in.nextInt());
			U3.add(in.nextInt());
		}
		in.close();
		Edge[] edges = new Edge[2*nEdges];
		
		int j=0;
		for(int i=0; i<(2*nEdges); i =i+2)
		{     		
			edges[i] = new Edge(U1.get(j), U2.get(j), U3.get(j));			
			edges[i+1] = new Edge(U2.get(j),U1.get(j),U3.get(j));
                 j=j+1;
		}
		
		return edges;
	}

//returning the no of nodes present */

	public static int getNumberOfNodes(Edge[] edges)
	{
		return n;
	}
	
	
        
	public static int[][] simpleScheme(Vertex[] v)
	{
		int n = v.length;
		int[][] edgeCost = new int[n][n];
		for(int i=0; i<n; i++)
			for(int j=0; j<n; j++)
				edgeCost[i][j] = Integer.MAX_VALUE;
		
		int Index = 0;
		int Dist;
		Set<Vertex> V = new HashSet<Vertex>();	/

/sets of nodes, the shortest paths to those nodes are not yet found
		
				
		
		for(int source=0; source<n; source++)
		{
			edgeCost[source][source] = 0;
			
			V.clear();
			
			
			
			for(Vertex vertex : v)
			{
				V.add(vertex);
			}
			
			
			while(!V.isEmpty())		

//find shortest path to n nodes

			{
				Dist = Integer.MAX_VALUE;
				
				for(Vertex gnode : V)	

//find the node that has the shortest distance among all undetermined nodes

				{
					if(Dist > edgeCost[source][gnode.getIndex()])
					{
						Dist = edgeCost[source][gnode.getIndex()];
						Index = gnode.getIndex();
					}
				}
				
				V.remove(v[Index]);
				
				
//the distance between the source and a neighbor of the new added node may be reduced

				for(neighbors an : v[Index].getNeighbors())
				{
					if(edgeCost[source][an.getIndex()] > edgeCost[source][Index] + an.getDistance())
						edgeCost[source][an.getIndex()] = (int) (edgeCost[source][Index] + an.getDistance());
					
				}	
			
			
			
			}
			
			
		}
		return edgeCost;
	}
	
	
	public static int[][] fheapScheme(Vertex[] g)
	{ 

/* generating the fibonacci heap*/

		FibonacciHeap fheap = new FibonacciHeap();
		int n = g.length;
		int[][] edgeCost = new int[n][n];
		for(int i=0; i<n; i++)
			for(int j=0; j<n; j++)
				edgeCost[i][j] = Integer.MAX_VALUE;
		int Index = 0;
		Set<Vertex> V = new HashSet<Vertex>();
	
		for(int source=0; source<n; source++)
		{
			edgeCost[source][source] = 0;
			V.clear();
			for(Vertex gnode : g)
			{
				V.add(gnode);
			}
/*inserting the source into the fibonacci heap*/

			fheap.insert(source, 0);
			
			V.remove(g[source]);
			
			while(fheap.getMin() != null)
			{
				Index = fheap.extractMin().getIndex();
							
//distance to the adjacent nodes of the recently removed node may be reduced
				for(neighbors an : g[Index].getNeighbors())
				{
					if(edgeCost[source][an.getIndex()] > edgeCost[source][Index] + an.getDistance())
					{
						edgeCost[source][an.getIndex()] = (int) (edgeCost[source][Index] + an.getDistance());

//if the node has not been inserted to the heap, then insert it to the heap

						if(V.contains(g[an.getIndex()]))
						{
							fheap.insert(an.getIndex(), edgeCost[source][an.getIndex()]);
							V.remove(g[an.getIndex()]);
						}
						
									
//if the node is already in the heap, find the corresponding heap node, and decrease the distance


						else
							fheap.decreaseKey(fheap.searchNode(an.getIndex(), fheap.getMin()), edgeCost[source][an.getIndex()]);
					}
				}	
			}
		}
		return edgeCost;
	}

	
	
	
	public static Vertex[] generateGraph(Edge[] edges, int n)
	{
		Vertex[] vertices = new Vertex[n];
		for(int i=0; i<n; i++)
		{
			vertices[i] = new Vertex(i);
		}
		
		for(Edge e: edges)
		{
			vertices[e.getVertex1()].insertNeighbor(e.getVertex2(), e.getDistance());
		}
			
		return vertices;
	}
	
	public static boolean isConnected(Vertex[] v)
	{
		int[][] edgeCost = simpleScheme(v);
		
//execute simple scheme, if the distance between any pair of nodes is 

//larger than the possible maximum distance if they are connected, then

//the graph is not connected
		
		for(int i=0; i<v.length; i++)
			for(int j=0; j<v.length; j++)
			{
				if(edgeCost[i][i] > 1000 * v.length + 1)
					return false;
			}
		return true;
	}

}

/*  vertex here refers to the nodes in the graph */

class Vertex 
{
	private int index;
	private ArrayList<neighbors> neighbors;		

//list of adjacent nodes
	
	public Vertex(int index)
	{
		this.index = index;
		neighbors = new ArrayList<neighbors>();
	}
	
	
	public void insertNeighbor(int i, int d)
	{
		neighbors.add(new neighbors(i, d));
	}
	
	
	public ArrayList<neighbors> getNeighbors()
	{
		return neighbors;
	}
	
	public int getIndex()
	{
		return index;
	}
	
}

class Edge 
{
  private int vertex1;
	private int vertex2;
	private int distance;
	
	
	public Edge(int start, int end, int distance)
	{
		this.vertex1 = start;
		this.vertex2 = end;
		this.distance = distance;
	}
	
	public int getVertex1()
	{
		return vertex1;
	}
	
	public int getVertex2()
	{
		return vertex2;
	}
	
	public int getDistance()
	{
		return distance;
	}
	
	public boolean equals(Edge e)
	{
		if(e != null && vertex1 == e.getVertex1() && vertex2 == e.getVertex2())
			return true;
		else return false;
	}


}


/* adjacency class : returning and keeping the track of next neighbours */

class neighbors 
{
  private int index;
	private double distance;
	
/* Constructor :     i node index : d distance to that adjacent node  */

	public neighbors(int i, double d)
	{
		index = i;
		distance = d;
	}
	
	public int getIndex()
	{
		return index;
	}

	public double getDistance()
	{
		return distance;
	}
}


class FibonacciHeap 
{
  private Fnode min;
	private int size;
	
	
 	public static class Fnode
	{
		private int index;
		private int degree;
		private double distance;
		private Fnode parent;
		private Fnode child;
		private Fnode leftSibling;
		private Fnode rightSibling;
		private boolean childCut;


/* the childcut field is set true iif the node has lost a child since it became child of its current parent and set to false by remove min, which is the only operation that makes one node a child of another */
		
		public Fnode(int index, double distance)
		{
			this.index = index;
			this.distance = distance;
			this.degree = 0;
			this.childCut = false;
		}
		
		public void setIndex(int i)
		{
			index = i;
		}
		
		public int getIndex()
		{
			return index;
		}
		
		public void setDistance(double d)
		{
			distance = d;
		}
		
		public double getDistcance()
		{
			return distance;
		}
		
		public void setParent(Fnode p)
		{
			parent = p;
		}
		
		public Fnode getParent()
		{
			return parent;
		}
		
		public void setChild(Fnode c)
		{
			child = c;
		}
		
		public Fnode getChild()
		{
			return child;
		}
		
		public void setLeftSibling(Fnode ls)
		{
			leftSibling = ls;
		}
		
		public Fnode getLeftSibling()
		{
			return leftSibling;
		}
		
		public void setRightSibling(Fnode rs)
		{
			rightSibling = rs;
		}
		
		public Fnode getRightSibling()
		{
			return rightSibling;
		}
		
		public void increaseDgree()
		{
			degree ++;
		}
		
		public void decraseDegree()
		{
			degree --;
		}
		
		public int getDegree()
		{
			return degree;
		}
		
		public void setChildCut(boolean b)
		{
			childCut = b;
		}
		
		public boolean getChildCut()
		{
			return childCut;
		}
		
	}

 /* FibonacciHeap() : constructor create an empty Fibonacci heap*/


	public FibonacciHeap()
	{
		min = null; //generating empty heap : initializing the minimum as null
		size = 0;
	}


/* generating the fibonacci heap and setting the left and right siblings */	
	
	public FibonacciHeap(Fnode x)
	{
		min = x;
		x.setChild(null);
		x.setParent(null);
		x.setLeftSibling(x);
		x.setRightSibling(x);
		size = 1;
	}
	
/* returning the minimun from the fibonacci heap */


	public Fnode getMin()
	{
		return min;
	}
	
/* returning the size of the heap */	


	public int getSize()
	{
		return size;
	}
	


/*the function link() is removing y from root list, and make y a child of x, y. it can be noted that distance must be larger than x.distance, also y a node removed from the root list, becomes x' child    x a node in the root list, becomes y' parent*/


	public void link(Fnode y, Fnode x) 		
	{											
		Fnode l = y.getLeftSibling();
		Fnode r = y.getRightSibling();
		
		l.setRightSibling(r);
		r.setLeftSibling(l);
		y.setParent(x);
		y.setChildCut(false);
		
		if(x.getDegree() == 0)				

//if y becomes the only child
		{
			y.setLeftSibling(y);
			y.setRightSibling(y);
		}
		
		else								

// add y to the child list of x
		{
			y.setRightSibling(x.getChild());
			y.setLeftSibling(x.getChild().getLeftSibling());
			x.getChild().getLeftSibling().setRightSibling(y);
			x.getChild().setLeftSibling(y);
		}
			
		x.setChild(y);
		x.increaseDgree();
			
	}
	

 //Union with another Fibonacci heap h
	 
	public void union(FibonacciHeap h)
	{
		Fnode min2 = h.getMin();
		size += h.getSize();
		
		//if(min2 == null) do nothing
		
		if(min == null && min2 !=null)
		{
			min = min2;
		}

		if(min != null && min2 != null)		

//add h's root list to this heap's root list, adjust min if necessary
		{
			Fnode l = min.getLeftSibling();
			
			min2.getLeftSibling().setRightSibling(min);
			l.setRightSibling(min2);
			min.setLeftSibling(min2.getLeftSibling());
			min2.setLeftSibling(l);
			
			if(min.getDistcance() > min2.getDistcance())
				min = min2;
		}
	}
	
							

// inserting the node into the heap


	public void insert(int index, double dist)
	{
		Fnode fnode = new Fnode(index, dist);
		FibonacciHeap h = new FibonacciHeap(fnode);
		this.union(h);
	}
	
	
	public void consolidate()
	{
								

//possible maximum degree for a F-heap with "size" nodes


		int bucketSize = (int) (Math.log(size) / Math.log(2)) + 2; 	
								
//degreeMap[i] means the node whose degree is i


		Fnode[] degreeMap = new Fnode[bucketSize];					
		for(int i=0; i<bucketSize; i++)
			degreeMap[i] = null;
		
		int d;
		Fnode x = min;
		Fnode start = min;
		Fnode y, next;
		do
		{
			d = x.getDegree();
			next = x.getRightSibling();
			while(degreeMap[d] != null)
			{
				y = degreeMap[d];
				if(x.getDistcance() > y.getDistcance())		


//exchange x, y if necessary (refer to link())


				{
					Fnode temp = y;
					y = x;
					x = temp;
				}
				if(y == start)								

//maintain the ending mark for the loop


					start = start.getRightSibling();
				
				if(y == next)								

//maintain the next pointer for x (should be in the root list)


					next = next.getRightSibling();
				this.link(y, x);
				degreeMap[d] = null;
				d ++;
			}
			degreeMap[d] = x;
			x = next;
		}while(x != start);
		
		
		min = null;
		for(int i=0; i<bucketSize; i++)	//adjust min 
		{
			if(min == null || (degreeMap[i] != null && min.getDistcance() > degreeMap[i].getDistcance()))
				min = degreeMap[i];
		}
	}
	
	
/* sets the childcut filed to false, and is the only operation that makes one node a child of another */

	public Fnode extractMin()
	{
		Fnode z = min;
		if(z != null)
		{
			Fnode x = z.getChild();		

// set children's parent to null

			if(x != null)
			{
				x.setParent(null);
				for(x = x.getRightSibling(); x != z.getChild(); x = x.getRightSibling())
					x.setParent(null);
			
			
				Fnode l = z.getLeftSibling();
				x = z.getChild();
			
				l.setRightSibling(x);	

//add z's children to root list

				x.getLeftSibling().setRightSibling(z);
				z.setLeftSibling(x.getLeftSibling());
				x.setLeftSibling(l);
			}
			
			
			z.getLeftSibling().setRightSibling(z.getRightSibling());	

//remove z from root list

			z.getRightSibling().setLeftSibling(z.getLeftSibling());
			
					
			
//this means that z is the only root even after it children are added to root list (no children)

			if(z == z.getRightSibling())		
				min = null;
				
			else
			{
				min = z.getRightSibling();
				this.consolidate();		

// adjustment of min is done in consolidate()

			}
			size--;
		}
		return z;
	}

/* decreaseKey() function removes subtree rooted at theNode from its doubly linked sibling list if theNode is not a root and new key < parent key */
	
	public void decreaseKey(Fnode x, double k)
	{
		if(k > x.getDistcance())
		{
			System.out.println("Error! New key is greater than current key.");
			System.exit(1);
		}
		
		x.setDistance(k);
		Fnode y = x.getParent();
		if(y != null && x.getDistcance() < y.getDistcance())
		{
			cut(x, y);			

//move x to root list

			cascadingCut(y);	

//cascading cut if necessary

		}
		
		if(x.getDistcance() < min.getDistcance())
			min = x;			

//adjust min if necessary

	}
	
	
	public void cut(Fnode x, Fnode y)  					

//y is x's parent

	{
		Fnode z = y.getChild();							

//remove x from the child list of y

		if (z == x && x.getRightSibling() == x)			//x is y's only child
			y.setChild(null);
		
		if (z == x && x.getRightSibling() != x)			//x is not the only child, but the pointer points to x
		{
			Fnode l = x.getLeftSibling();
			Fnode r = x.getRightSibling();
			l.setRightSibling(r);
			r.setLeftSibling(l);
			y.setChild(r);
		}
		
		if(z != x)				//x is not the only child, and the pointer does not point to x
		{
			Fnode l = x.getLeftSibling();
			Fnode r = x.getRightSibling();
			
			l.setRightSibling(r);
			r.setLeftSibling(l);
		}
		
		y.decraseDegree();
		
		x.setRightSibling(min);						

//add x to root list, as y exists, min is not null

		min.getLeftSibling().setRightSibling(x);
		x.setLeftSibling(min.getLeftSibling());
		min.setLeftSibling(x);
		
		x.setParent(null);
		x.setChildCut(false);
	}
	
	
	public void cascadingCut(Fnode y)
	{
		Fnode z = y.getParent();
		if(z != null)
		{
			if(y.getChildCut() == false)
				y.setChildCut(true);
			
			else
			{
				this.cut(y, z);
				this.cascadingCut(z);
			}
		}
	}
	
	
	public void delete(Fnode x)
	{
		this.decreaseKey(x, Double.NEGATIVE_INFINITY);
		this.extractMin();
	}
	
	/
	public Fnode searchNode(int index, Fnode x) 
	{
		Fnode wantedNode = null;
		if(x != null)
		{
			Fnode y  = x;
			do
			{
				if(y.getIndex() == index)
					return y;
				if((wantedNode = searchNode(index, y.getChild())) != null)
					return wantedNode;
				y = y.getRightSibling();
				
			}while (y != x);
		}
		return wantedNode;	
	}
	
	
	public void print(Fnode x)
	{
		Fnode y = x;
		if(x != null)
		{
			do
			{
				System.out.print(y.getIndex() + "\t");
				print(y.getChild());
				
				System.out.println();
				
				y = y.getRightSibling();
			}
			while(y != x);
		}
	}
	
}

