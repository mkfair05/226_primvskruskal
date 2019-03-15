/* PrimVsKruskal.java
   CSC 226 - Spring 2019
   Assignment 2 - Prim MST versus Kruskal MST Template
   
   The file includes the "import edu.princeton.cs.algs4.*;" so that yo can use
   any of the code in the algs4.jar file. You should be able to compile your program
   with the command
   
	javac -cp .;algs4.jar PrimVsKruskal.java
	
   To conveniently test the algorithm with a large input, create a text file
   containing a test graphs (in the format described below) and run
   the program with
   
	java -cp .;algs4.jar PrimVsKruskal file.txt
	
   where file.txt is replaced by the name of the text file.
   
   The input consists of a graph (as an adjacency matrix) in the following format:
   
    <number of vertices>
	<adjacency matrix row 1>
	...
	<adjacency matrix row n>
	
   Entry G[i][j] >= 0.0 of the adjacency matrix gives the weight (as type double) of the edge from 
   vertex i to vertex j (if G[i][j] is 0.0, then the edge does not exist).
   Note that since the graph is undirected, it is assumed that G[i][j]
   is always equal to G[j][i].


   R. Little - 03/07/2019
*/

 import edu.princeton.cs.algs4.*;
 import java.util.Scanner;

//  import org.graalvm.compiler.lir.EdgeMoveOptimizer;

 import java.io.File;

//Do not change the name of the PrimVsKruskal class
public class PrimVsKruskal{

	static Queue<Edge> primMST;  //mark edges in Prim MST
	static boolean[] marked;     //mark visited edges in prim's
	static double[] distTo;		 // weight of shortest edge
	static Edge[] edgeTo;	 	 //shortest edge from tree vertex to non tree vertex
	static IndexMinPQ<Double> primPQ;   //priority queue for prim implementation
	//-----------------------------------------------------------------------
	static Queue<Edge> krusMST;  //mark edges in Kruskal MST
	static IndexMinPQ<Edge> krusPQ;   //priority queue for kruskal implementation
	static UF uf;				  //Union find data structure
	static double weight;		 //total weight of kruskal MST
	//-----------------------------------------------------------------------
	//gets the adjacency list from the matrix
	public static double[] GetAdjList(double[][] G, int index){
		return G[index];
	}

	//prints adjacency list
	static void PrintAdjList(double[] list){
		for (int i = 0; i < list.length; i++){
			System.out.println(list[i]);
		}
	}

	static int GetLength(double[][] G){
		return G.length;
	}

	static void InitPrim(double[] adjList, double[][] G){
		int numVertices = GetLength(G);
		
		primPQ = new IndexMinPQ<Double>(numVertices);
		edgeTo = new Edge[numVertices];
		distTo = new double[numVertices];
		marked = new boolean[numVertices];

		for (int i = 0; i < numVertices; i++){
			distTo[i] = Double.POSITIVE_INFINITY;
		}
		for (int i = 0; i < numVertices; i++){
			if (!marked[i]){
				RunPrims(adjList, i, G);
			}
		}
	}

	static void RunPrims(double[] adjList, int vertex, double[][] G){
		distTo[vertex] = 0.0;
		primPQ.insert(vertex, distTo[vertex]);
		while(!primPQ.isEmpty()){
			int m = primPQ.delMin();
			ScanVertex(adjList, m, G); //scan all the edges adjacent to m
		}
	}

	static void ScanVertex(double[] adjList, int vertex, double[][] G){
		marked[vertex] = true;
		int numVertices = GetLength(G);

		for (int i = 0; i < numVertices; i++){
			if (adjList[i] > 0){
				//there exists an edge between vertices
				int w = i;
				if (marked[w]){
					//vertex and w is obsolete edge
					continue;
				}
				if (adjList[i] < distTo[w]){
					distTo[w] = adjList[i];
					edgeTo[w] = new Edge(vertex, w, adjList[i]);
					if (primPQ.contains(w)){
						primPQ.decreaseKey(w, distTo[w]);
					} else {
						primPQ.insert(w, distTo[w]);
					}
				}
			}
		}
	}

	static void InitKruskal(double[] adjList, int vertex, double[][] G){
		krusPQ = new MinPQ<Edge>();
		int numVertices = GetLength(G);
		for (int i = 0; i < adjList.length; i++){ //get all edges around A
			if (adjList[i] > 0){
				//there exists an edge between vertices index and i
				int w = i;  
				double weight = adjList[i];
				Edge e = new Edge(vertex, w, weight);
				krusPQ.insert(e); //insert all edges into priority queue
			}
		}

		UF uf = new UF(numVertices);
		while (!krusPQ.isEmpty() && krusMST.size() < numVertices){
			Edge e = krusPQ.delMin();
			int v = e.either();
			int w = e.other(v);
			if (!uf.connected(v,w)){
				//if v-w doesn't create a cycle
				uf.union(v,w);
				krusMST.enqueue(e); //add edge to MST
				weight +=e.weight();
			}
		}
	}

	/* PrimVsKruskal(G)
		Given an adjacency matrix for connected graph G, with no self-loops or parallel edges,
		determine if the minimum spanning tree of G found by Prim's algorithm is equal to 
		the minimum spanning tree of G found by Kruskal's algorithm.
		
		If G[i][j] == 0.0, there is no edge between vertex i and vertex j
		If G[i][j] > 0.0, there is an edge between vertices i and j, and the
		value of G[i][j] gives the weight of the edge.
		No entries of G will be negative.
	*/
	static boolean PrimVsKruskal(double[][] G){
		/* Build the MST by Prim's and the MST by Kruskal's */
		/* (You may add extra methods if necessary) */
		
		/* ... Your code here ... */
		
		int n = GetLength(G);
		double[] adjList;
		for (int i = 0; i < n; i++){
			adjList = GetAdjList(G, i);
			InitKruskal(adjList, i, G);
			InitPrim(adjList, G);
		}










		
		/* Determine if the MST by Prim equals the MST by Kruskal */
		boolean pvk = true;
		/* ... Your code here ... */

		return pvk;	
	}
		
	/* main()
	   Contains code to test the PrimVsKruskal function. You may modify the
	   testing code if needed, but nothing in this function will be considered
	   during marking, and the testing process used for marking will not
	   execute any of the code below. 
	*/
   public static void main(String[] args) {
		Scanner s;
		if (args.length > 0){
			try{
				s = new Scanner(new File(args[0]));
			} catch(java.io.FileNotFoundException e){
				System.out.printf("Unable to open %s\n",args[0]);
				return;
			}
			System.out.printf("Reading input values from %s.\n",args[0]);
		}else{
			s = new Scanner(System.in);
			System.out.printf("Reading input values from stdin.\n");
		}
		
		int n = s.nextInt();
		double[][] G = new double[n][n];
		int valuesRead = 0;
		for (int i = 0; i < n && s.hasNextDouble(); i++){
			for (int j = 0; j < n && s.hasNextDouble(); j++){
				G[i][j] = s.nextDouble();
				if (i == j && G[i][j] != 0.0) {
					System.out.printf("Adjacency matrix contains self-loops.\n");
					return;
				}
				if (G[i][j] < 0.0) {
					System.out.printf("Adjacency matrix contains negative values.\n");
					return;
				}
				if (j < i && G[i][j] != G[j][i]) {
					System.out.printf("Adjacency matrix is not symmetric.\n");
					return;
				}
				valuesRead++;
			}
		}
		
		if (valuesRead < n*n){
			System.out.printf("Adjacency matrix for the graph contains too few values.\n");
			return;
		}	
        boolean pvk = PrimVsKruskal(G);
        System.out.printf("Does Prim MST = Kruskal MST? %b\n", pvk);
    }
}
