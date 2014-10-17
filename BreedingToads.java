import java.util.*;
import java.io.*;
public class BreedingToads{

   private static final int POINTS = 100;

   static double Xa,Xb,Ya,Yb = 0;

   public static void main(String[]args){
     try{
      Scanner s  = new Scanner(new File("./Toads.txt"));    
      double[][] pts = readData(s);
      int[] subSet = getMinKthNearestN(pts,13);
      double[][] bounds = getBounds(pts,subSet);
      double[] centroid = getCentroid(pts,subSet);
      double [] step = {1,1};
      double r = scan(pts,bounds,step);
      System.out.println("largest Radius: " + r);
     }catch(Exception e){
         System.err.println("input.txt aint there");
         System.err.println(e);
      }
   }                      

   public static int[] findClosestTo(double[][] pts,double[] centroid,int[] candidates){
      double closestSoFar = Math.pow(centroid[0]-pts[candidates[0]][0],2)+Math.pow(centroid[1]-pts[candidates[0]][1],2);
      int[] pos = new int[1];
      for(int i = 1;i<candidates.length;i++){
         double distance = Math.pow(centroid[0]-pts[candidates[i]][0],2)+Math.pow(centroid[1]-pts[candidates[i]][1],2);
         if(distance < closestSoFar){
            pos[0] = i;
         }
      }
      return pos;
   }

   public static double[] getCentroid(double[][] pts, int[] subSet){
                                                                                                              
      double minX = Double.MAX_VALUE, maxX = 0, minY = Double.MAX_VALUE, MaxY = 0;
                                                                                                            
      //iterate over subSet        
      double x_sum = 0.0;
      double y_sum = 0.0;

      //add our sums in pts defind by the indexes at subSet
      for(int i = 0; i < subSet.length; i++){                          
         x_sum += pts[subSet[i]][0];
         y_sum += pts[subSet[i]][1];           
      }
    //get the average of these which is our centroid
      x_sum /= subSet.length;
      y_sum /=subSet.length;
 
      //return the two values as a single entry in an array!
      double[] return_array = new double[2];        
      return_array[0] = x_sum;
      return_array[1] = y_sum;
      return return_array;
   }


   public static double[][] getBounds(double[][] pts){
      int[] subSet = new int[pts.length];
      for(int i = 0;i<subSet.length;i++){
         subSet[i] = i;
      }
      return getBounds(pts,subSet);
   }
                                                                                                                             
   /**
    * Returns the bounds of a Subset of pts defined by subSet.
    *
    * @param pts Array of points in 2-d Space where pts[i][j] refers
    * to the ith point and j is the dimension {j=0:x,j=1:y}
    * @param subSet Indexs of pts of which the subset is made up of.
    * @returns A set of points which describe the bounding box
    * {{x1,y1},{x2,y2}} where 1 is one corner and 2 is the diagonal
    * opposite of the bounding box of the Sub set.
    */
   public static double[][] getBounds(double[][] pts, int[] subSet){
     double minX = Double.MAX_VALUE;
      double maxX = 0;
      double minY = Double.MAX_VALUE;
      double maxY = 0;
      //iterate over subSet
      for(int i = 0; i < subSet.length; i++){                           
         if(pts[subSet[i]][0] < minX){
            minX = pts[subSet[i]][0];
         }
         if(pts[subSet[i]][0] > maxX){
            maxX = pts[subSet[i]][0];
         }       
         if(pts[subSet[i]][1] > maxY){
            maxY = pts[subSet[i]][1];
         }
         if(pts[subSet[i]][1] < minY){
            minY = pts[subSet[i]][1];
         }
      }           
      double[][] return_array = new double[2][2];


      return_array[0][0] = minX;
      return_array[0][1] = minY;
      return_array[1][0] = maxX;
      return_array[1][1] = maxY;

      return return_array;
   }
          
   public static double[] getStep(double[][]pts, int[] subSet){
      double[] dist = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};   
      for(int i = 0; i<subSet.length;i++){
         for(int j = 0; j<subSet.length;j++){
            double distX =  pts[subSet[i]][0] - pts[subSet[j]][0];
            double distY =  pts[subSet[i]][1] - pts[subSet[j]][1];
            if(dist[0]>distX){
               dist[0] = distX;
            }
            if(dist[1] >distY){
               dist[1] = distY;
            }
         }
      }
      // We want to step at less the minimum distance between points.
      dist[0] = Math.abs(dist[0]/2);
      dist[1] = Math.abs(dist[0]/2);
      if(dist[0] == 0){
         dist[0] = 0.0001;
      }
      if(dist[1] == 0){
         dist[1] = 0.0001;
      }
      return dist;
   }
                                                                                                                                                     
   /**
    * Parses from the Scanner parameter co-ordinates and returns them in a double[][].
    *
    * @param scan Scanner which contains a Header line "telephone sites" then pairs of doubles.
    *
    * @return double[i][j] where i runs over the points and j =0 is x
    * and j=1 is y. Assuming x values preceed y values in scan.
    *
    */
   public static double[][] readData(Scanner scan){
                              
      double maxX = 0;
      double maxY = 0;
      ArrayList<String> coords = new ArrayList<String>();
      while(scan.hasNextLine()){
            coords.add(scan.nextLine());
      }
      Scanner stringParser;
      double[][] pts = new double[coords.size()][2];
      for(int i = 0; i < coords.size(); i++){
         String line = coords.get(i);  
         stringParser = new Scanner(line);
         pts[i][0] = stringParser.nextDouble(); 
         if(pts[i][0] > maxX){
            maxX = pts[i][0];
         }
         pts[i][1] = stringParser.nextDouble();          
         if(pts[i][1] > maxY){
            maxY = pts[i][1];
         }   
      }
      return pts;
   }
                                       
   public static void eprintNearestNeighbours(double[][][] nearestN){
      for(int i = 0; i< nearestN.length; i++){
         for(int j = 0; j < nearestN[0].length;j++){
            System.err.print(j + " Nearest Neighbour to " + i + " is "  + nearestN[i][j][0] + " distance "
                             + nearestN[i][j][1] + "\n");
         }
         System.err.print("\n");
      }
   }
                                                        
   /**
    * Finds the set of points that are the closest together.
    *
    * @param pts Set of points to search.
    * @param k Number of points required in the subset.
    *
    * @returns An array of the index's in pts which make up the subset.
    */
   public static int[] getMinKthNearestN(double[][] pts, int k){
      if(k>= pts.length){
         int[] indexs = new int[k];
         for(int i = 0;i<indexs.length;i++){
            indexs[i] = i;
         }
         return indexs;
      }

      double[][] distanceTo = new double[pts.length][pts.length];
      int[][] indexs = new int[pts.length][pts.length];
      for(int i = 0; i<pts.length;i++){
         for(int j = 0; j<pts.length;j++){
            indexs[i][j] = j;
            distanceTo[i][j] = Math.pow(pts[i][0] - pts[j][0],2)+Math.pow(pts[i][1] - pts[j][1],2);
         }
      }


      for(int i = 0; i<pts.length;i++){
         iSort(indexs[i],distanceTo[i]);
      }

      double minKthNearestN = Double.POSITIVE_INFINITY;
      int pos = -1;
      for(int i = 0; i<pts.length;i++){
         if(distanceTo[i][k] < minKthNearestN){
            minKthNearestN = distanceTo[i][k];
            pos = i;
         }
      }
      int[] space = new int[k+1];
      for(int i = 0; i<=k; i++){
         space[i] = indexs[pos][i];
      }

      return space;
   }

   /**
    * iSort
    *
    * Use insertion sort to sort the arrays keys and values based on
    * values.
    * @param keys Auxillary array which is sorted based on values in values
    * @param values Array of values which are sorted.
    */
   public static void iSort(int[]keys, double[] values){
      for(int j = 0;j<keys.length;j++){
         int k = keys[j];
         double v = values[j];
         int i = j - 1;
         while((i>=0)&&(values[i]>v)){
            keys[i+1] = keys[i];
            values[i+1] = values[i];
            i--;
         }
         keys[i+1] = k;
         values[i+1] = v;
      }
   }

   public static int[] getPtsIn(double[][] pts,double[][] boundingBox){
      ArrayList<Integer>subSet = new ArrayList<Integer>();
      for(int i = 0;i<pts.length;i++){
         subSet.add(new Integer(i));

      }
      int[] sSet = new int[subSet.size()];
      for(int i=0;i<sSet.length;i++){
         sSet[i] = subSet.get(i).intValue();
      }
      return sSet;
   }   

   /**Scans across search space finding the largest bounding circle
    * of no more than 12 toads.
    *
    * @param pts Points testing for
    * @param searchSpace The space to search over. {{x1, y1},{x2,y2}}.
    * @param step The step size of the search step space. {xStep, yStep}
    */
   public static double scan(double[][] pts, double[][] searchSpace, double[] step){   
      double r = Double.POSITIVE_INFINITY;        
      if(pts.length<13){
         return r;
      }       
      double deltaX = (searchSpace[1][0]-searchSpace[0][0]);
      double deltaY = (searchSpace[1][1]-searchSpace[0][1]);  
      int []subSet;
      do{
         double[][] boundingBox = new double[2][2];
         boundingBox[0][0] = searchSpace[0][0]-deltaX;
         boundingBox[1][0] = searchSpace[1][0]+deltaX;
         boundingBox[0][1] = searchSpace[0][1]-deltaY;
         boundingBox[1][1] = searchSpace[1][1]+deltaY;
         subSet = getPtsIn(pts,boundingBox);
         deltaX*=2;
         deltaY*=2;
      }while(subSet.length<13);

      final int MAX_ENCLOSED = 12;

      for(double x = searchSpace[0][0];x<=searchSpace[1][0];x+=step[0]){                
         for(double y = searchSpace[0][1];y<=searchSpace[1][1];y+=step[1]){          
            double[] distances = new double[pts.length];
            for(int i = 0; i < subSet.length; i++){
               distances[i] = Math.pow(pts[subSet[i]][0] - x,2) + 
                  Math.pow(pts[subSet[i]][1]-y,2);
            }
            Arrays.sort(distances);
            if(r>distances[12]){
               r = distances[12];
            }       
         }
      }
      System.out.println();
      return Math.sqrt(r);
   }   
}



