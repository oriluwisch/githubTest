import java.util.*;

public class Graph {
    public Vert [] verts;
    public Set<Edge> edges;C:\Users\ori\IdeaProjects\WI-research\src\Graph.java
    public int numOfVertices;
    public int numOfEdges;

    public Graph(){}
    public Graph (Vert []_verts, Set<Edge> _edges){
        verts = _verts;
        edges = _edges;
        numOfVertices = verts.length;
        numOfEdges = edges.size();
    }

    public void clearEdges(){
        edges.clear();
        numOfEdges = 0;
    }

    public static double[] linearEq(Vert v1, Vert v2){
        double [] coef = new double[2];
        double x1 = v1.coorX;
        double y1 = v1.coorY;
        double x2 = v2.coorX;
        double y2 = v2.coorY;
        coef[0] = Math.abs((y2-y1)/(x2-x1));
        coef[1] = y1 - coef[0]*x1;
        return coef;
    }


    public void print2D(double[][] mat){
        for(int i = 0; i < mat[0].length; i++){
            System.out.println();
            for(int j = 0; j < mat[0].length; j++){
                System.out.print((float) Math.round(mat[i][j] * 100) / 100 + " ");
            }
        }
    }
    public double[][] getDistMat(){
        double[][] distMat = new double[numOfVertices][numOfVertices];
        for(int i = 0; i < numOfVertices; i++){
            for(int j = 0; j < numOfVertices; j++){
                if(i != j) distMat[i][j] = 10000;
                else distMat[i][j] = 0;
            }
        }
//        print2D(distMat);
//        System.out.println();
        for(Edge edge: edges){
            distMat[edge.getStartVert().index][edge.getTargetVert().index] = edge.weight;
            distMat[edge.getTargetVert().index][edge.getStartVert().index] = edge.weight;
        }
//        print2D(distMat);
//        System.out.println();

        for(int k = 0; k < numOfVertices; k++){
            for(int j = 0; j < numOfVertices; j++){
                for(int i = 0; i < numOfVertices; i++) {
                    if (distMat[i][k] + distMat[k][j] < distMat[i][j])
                        distMat[i][j] = distMat[i][k] + distMat[k][j];
                }
            }
        }
//        print2D(distMat);
//        System.out.println();
        return distMat;
    }

    public double getWI(){
        double [][] distMat = getDistMat();
        double WI = 0.0;
        for(int i = 0; i < distMat[0].length; i++){
            for(int j = 0; j < distMat[0].length; j++){
                WI += distMat[i][j];
            }
        }
        return 0.5 * WI;
    }

    public void printEdge(Edge e){
        System.out.print("("+e.startVert.getName() + "," + e.targetVert.getName() + ") ");
    }

    public void printEdges(Set<Edge> edges){
        int i = 0;
        for(Edge e: edges) {
            printEdge(e);
            if(i%5==1) System.out.println();
            i++;
        }
    }

    public static boolean closer(Vert v1, Vert v2, Vert v3){
        if (v1.euclideanDistance(v3) < v2.euclideanDistance(v3)) return true;
        return false;
    }

    public static void fwTest(){
        Vert[] edgeCaseVerticesArr = new Vert[300];
        for (int i = 0; i < edgeCaseVerticesArr.length; i++) {
            edgeCaseVerticesArr[i] = new Vert("v_" + String.valueOf(i), i);
        }
//        for (int i = 0; i < edgeCaseVerticesArr.length; i += 2) {
//            edgeCaseVerticesArr[i] = new Vert("v_" + String.valueOf(i), 0, 1.0 + 0.000001 * i, i);
//        }
//        for (int i = 1; i < edgeCaseVerticesArr.length; i += 2) {
//            edgeCaseVerticesArr[i] = new Vert("v_" + String.valueOf(i), 1, 1.0 + 0.000001 * i, i);
//        }

        double completeGraphWI = 0;
        double starGraphWI = 0;

        // Calculate the WI when the graph is star graph from the point that minimize sum of distance
        int index = Vert.findStarCenter(edgeCaseVerticesArr);
        Set<Edge> edges = new HashSet<Edge>();
        for (int i=0; i < edgeCaseVerticesArr.length; i++){
            if(i != index) edges.add(new Edge(edgeCaseVerticesArr[i], edgeCaseVerticesArr[index]));
        }

        Set<Edge> edges1 = new HashSet<Edge>();
        for(int i = 0; i < edgeCaseVerticesArr.length; i++){
            for(int j = i+1; j < edgeCaseVerticesArr.length; j++){
                edges1.add(new Edge(edgeCaseVerticesArr[i], edgeCaseVerticesArr[j]));
            }
        }


        Graph testGraph = new Graph(edgeCaseVerticesArr, edges);
        testGraph.numOfVertices = edgeCaseVerticesArr.length;
        Graph testGraph1 = new Graph(edgeCaseVerticesArr, edges1);
        testGraph1.numOfVertices = edgeCaseVerticesArr.length;

        //testGraph.printEdges(testGraph.edges);
        for (int i = 0; i < edgeCaseVerticesArr.length; i++) {
            starGraphWI += edgeCaseVerticesArr[index].euclideanDistance(edgeCaseVerticesArr[i]);
        }
        starGraphWI *= (edgeCaseVerticesArr.length - 1);

        // Calculate the WI when the graph is the complete graph
        for (int i = 0; i < edgeCaseVerticesArr.length; i++) {
            for (int j = i + 1; j < edgeCaseVerticesArr.length; j++) {
                completeGraphWI += edgeCaseVerticesArr[i].euclideanDistance(edgeCaseVerticesArr[j]);
            }
        }

        System.out.println("The WI of testGraph by the floyd warshall calculation is " + testGraph.getWI());
        System.out.println("The WI of testGraph by the floyd warshall calculation is " + testGraph1.getWI());
        System.out.println("Wiener Index calculation- the two clusters case (with " + edgeCaseVerticesArr.length + " points):");
        System.out.println("The Wiener Index when the graph is the star graph is: " + starGraphWI);
        System.out.println("The Wiener Index when the graph is the complete graph is: " + completeGraphWI);
        System.out.println("The star graph gives Wiener Index with approximation of: " + starGraphWI / completeGraphWI);

    }

    public static void twoStarsTest(int Vsize){
        //System.out.println("Wiener Index calculation- the two stars case (with " + Vsize + " points):");
        /*Init n vertices randomly*/
        Vert[] verts = new Vert[Vsize];
        for (int i = 0; i < verts.length; i++) {
            verts[i] = new Vert("v_" + String.valueOf(i), i);
        }
        //Option- 2-approximate case
//        for (int i = 0; i < verts.length; i+=2) {
//            verts[i] = new Vert("v_" + String.valueOf(i), 0, 1+0.000001*i ,i);
//        }
//        for (int i = 1; i < verts.length; i+=2) {
//            verts[i] = new Vert("v_" + String.valueOf(i), 1, 1+0.000001*i ,i);
//        }

        // Calculate the WI when the graph is the complete graph
        double completeGraphWI = 0;
        for (int i = 0; i < verts.length; i++) {
            for (int j = i+1; j < verts.length; j++) {
                completeGraphWI += verts[i].euclideanDistance(verts[j]);
            }
        }

        // Calculate the WI when the graph is star graph from the point that minimize sum of distance
        double starGraphWI = 0.0;
        int index = Vert.findStarCenter(verts);
        for (int i = 0; i < verts.length; i++) {
            starGraphWI += verts[index].euclideanDistance(verts[i]);
        }
        starGraphWI *= (verts.length - 1);
        if(starGraphWI > 2*completeGraphWI){
            System.out.println("Error");
            for (int i = 0; i < verts.length; i++) {
                for (int j = i+1; j < verts.length; j++) {
                    System.out.print(verts[i].euclideanDistance(verts[j]) + ", ");
                }
            }
            System.out.println();
            for (int i = 0; i < verts.length; i++) {
                System.out.print(verts[i].euclideanDistance(verts[index]) + ", ");
            }
            System.exit(1);
        }

        Set<Edge> edges = new HashSet<Edge>();
        Graph testGraph = new Graph(verts, edges);

        for (int i=0; i < verts.length; i++){
            if(i != index) edges.add(new Edge(verts[i], verts[index]));
        }


        double currentTwoStarsWI = Integer.MAX_VALUE;
        double minTwoStarsWI = Integer.MAX_VALUE;
        Set<Vert> rightVerts= new HashSet<>();
        Set<Vert> leftVerts= new HashSet<>();
        double [] starsValues = new double[verts.length];
        for(int a = 0; a < verts.length; a++){
            for(int b = a+1; b < verts.length; b++){

                //Clear the edges and the partition of the plane
                rightVerts.clear();
                leftVerts.clear();

                //split the vertices to two subsets according to their positions w.r.t the line between v_a,v_b
                rightVerts.add(verts[a]);
                rightVerts.add(verts[b]);
                double[] coef = linearEq(verts[a], verts[b]);
                for(int j = 0; j < verts.length; j++){
                    if(j != a && j != b){
                        if(coef[0]*verts[j].coorX  - verts[j].coorY + coef[1] >= 0){
                            rightVerts.add(verts[j]);
                        }
                        else leftVerts.add(verts[j]);
                    }
                }

                //Compute the star for all vertex in this partition
                for (Vert uRight : rightVerts) starsValues[uRight.index] = Vert.computeStar(rightVerts, uRight);
                for (Vert vLeft: leftVerts) starsValues[vLeft.index] = Vert.computeStar(leftVerts, vLeft);

                //check all 2-stars graph possible for this partition of the plane
                double s1 = leftVerts.size();
                double s2 = rightVerts.size();
                double sRight, sLeft = 0.0;
                for (Vert uRight : rightVerts) {
                    for (Vert vLeft : leftVerts) {
//                        testGraph.clearEdges();
//                        for (Vert u : rightVerts) {
//                            if (!u.equals(uRight)) testGraph.edges.add(new Edge(uRight, u));
//                        }
//                        for (Vert v : leftVerts) {
//                            if (!v.equals(vLeft)) testGraph.edges.add(new Edge(vLeft, v));
//                        }
//                        testGraph.edges.add(new Edge(uRight, vLeft));
//                        double expWI = testGraph.getWI();
                            sRight = starsValues[uRight.index];
                            sLeft = starsValues[vLeft.index];
                            double e = uRight.euclideanDistance(vLeft);
                            //TODO - Compute stars for every u,v for this partition
                            currentTwoStarsWI = (s1-1)*sLeft + (s2-1)*sRight + s2*sLeft + s1*sRight + s1*s2*e;
//                            System.out.println((int)currentTwoStarsWI == (int)expWI);
                            if (currentTwoStarsWI < minTwoStarsWI) {
                                minTwoStarsWI = currentTwoStarsWI;
                                //System.out.println("WI = " + minTwoStarsWI + " with partition of (" + leftVerts.size()+","+rightVerts.size()+")");
                                //testGraph.printEdges(edges);
                                //System.out.println();
                            }
                    }
                }
            }
        }

        System.out.println("The Wiener Index when the graph is the 1-star graph is ~ " + (int)starGraphWI);
        System.out.println("The Wiener Index when the graph is the complete graph is ~ " + (int)completeGraphWI);
        System.out.println("The best WI when the graph is 2-stars graph is ~ " + (int)minTwoStarsWI);
        System.out.println("The star graph gives Wiener Index with approximation of: " + starGraphWI / completeGraphWI);
        System.out.println("The 2-star graph gives Wiener Index with approximation of: " + minTwoStarsWI / completeGraphWI);
    }



    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        for (int i = 0;  i < 50; i++) {
            twoStarsTest(50);
            System.out.println();
        }
        long end = System.currentTimeMillis();
        long elapsedTime = end - start;
        System.out.println("Time for this program is: " + elapsedTime);

    }
}

































class Vert{

    private boolean visited;
    private String name;
    private java.util.List<Edge> List;
    private double dist = Double.MAX_VALUE;
    private Vert pr;
    public double coorX;
    public double coorY;
    int index;



    public Vert(String name, int _index) {
        this.name = name;
        this.List = new ArrayList<>();
        Random rand = new Random();
        this.coorX = 100 * rand.nextDouble();
        this.coorY = 100 * rand.nextDouble();
        this.index = _index;

    }

    public Vert(String name, double x, double y, int _index) {
        this.name = name;
        this.List = new ArrayList<>();
        this.coorX = x;
        this.coorY = y;
        this.index = _index;
    }

    public static int findStarCenter(Vert [] vertArr){
        int ind = 0;
        int minSumOfDist = Integer.MAX_VALUE;
        int currSumOfDist = 0;
        for (int i = 0; i < vertArr.length; i++){
            currSumOfDist = 0;
            for (int j = 0; j < vertArr.length; j++){
                currSumOfDist += vertArr[i].euclideanDistance(vertArr[j]);
            }
            if(currSumOfDist < minSumOfDist){
                ind = i;
                minSumOfDist = currSumOfDist;
            }
            //System.out.println("The sum of distances where " + vertArr[i].getName() + " is the center is "+ currSumOfDist);
        }
        //System.out.println();
        return ind;
    }

    public static Vert findStarCenter(HashSet<Vert> vertArr){
        Vert center = null;
        int minSumOfDist = Integer.MAX_VALUE;
        int currSumOfDist = 0;
        for (int i = 0; i < vertArr.size(); i++){
            currSumOfDist = 0;
            for(Vert u: vertArr){
                for(Vert v: vertArr){
                    currSumOfDist+= u.euclideanDistance(v);
                }
                if(currSumOfDist < minSumOfDist){
                    center = u;
                    minSumOfDist = currSumOfDist;
                }
            }
        }
        return center;
    }

    public static double computeStar(Set<Vert> vertices, Vert u){
        double starValue = 0.0;
        for (Vert v: vertices){
            starValue += u.euclideanDistance(v);
        }
        return starValue;
    }

    public double euclideanDistance(Vert other){
        return Math.sqrt(Math.pow(this.coorX-other.coorX, 2) + Math.pow(this.coorY-other.coorY, 2));
    }



    public void printVert(){
        System.out.print( "(" + this.coorX + ", " + this.coorY + "), ");
    }

    public List<Edge> getList() {
        return List;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }



    public void setList(List<Edge> List) {
        this.List = List;
    }

    public void addNeighbour(Edge edge) {
        this.List.add(edge);
    }

    public boolean Visited() {
        return visited;
    }

    public void setVisited(boolean visited) {
        this.visited = visited;
    }

    public Vert getPr() {
        return pr;
    }

    public void setPr(Vert pr) {
        this.pr = pr;
    }

    public double getDist() {
        return dist;
    }

    public void setDist(double dist) {
        this.dist = dist;
    }

}


class Edge {
    double weight;
    public Vert startVert;
    public Vert targetVert;

    public Edge(Vert startVert, Vert targetVert) {
        this.weight = startVert.euclideanDistance(targetVert);
        this.startVert = startVert;
        this.targetVert = targetVert;
    }

    public void printEdge(Edge e){
        System.out.println("("+this.startVert.getName() + "," + this.targetVert.getName() + ")");
    }


    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

    public Vert getStartVert() {
        return startVert;
    }

    public void setStartVert(Vert startVert) {
        this.startVert = startVert;
    }

    public Vert getTargetVert() {
        return targetVert;
    }

    public void setTargetVert(Vert targetVert) {
        this.targetVert = targetVert;
    }

}
