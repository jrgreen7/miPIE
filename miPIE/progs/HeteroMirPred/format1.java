
import java.io.*;
import java.util.*;

class format1
{

  public static void main(String args[])
  {
     try
     {
         BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("21_features.txt")));
         PrintWriter p = new PrintWriter("21.feature.formatted");

         String data = "";

         while((data = r1.readLine()) != null)
         {
            StringTokenizer st = new StringTokenizer(data);
            String s = "1 ";
 
            for(int i=1; i<=20; i++)
             s += i+":"+st.nextToken()+" ";
             
            p.println(s);
         }   
            
        r1.close();
        p.close();     
         
     
     }

     catch(Exception e)
     {
         e.printStackTrace();
     }
  }

}

                      
   


