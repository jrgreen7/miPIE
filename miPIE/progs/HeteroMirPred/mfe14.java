//calculates the MFEI1 and MFEI4.


import java.io.*;
import java.util.*;

class mfe14
{

  public static void main(String args[])
  {
     try
     {
         BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream(args[0])));
                  
         PrintWriter p = new PrintWriter(args[0]+".mfe14");

         String data = null;
        
           
         r1.readLine(); //to remove the header        
 
         while((data = r1.readLine()) != null)
         {
                StringTokenizer st = new StringTokenizer(data);
                String s = st.nextToken()+ "  "; //ID

                for(int i=1; i<=27; i++)
                  st.nextToken();
                                      
                double cg = Double.parseDouble(st.nextToken());

                for(int i=1; i<=17; i++)
                  st.nextToken(); 
                    
                double bp = Double.parseDouble(st.nextToken());
                st.nextToken();
                double mfe = Double.parseDouble(st.nextToken());
                double Nmfe = Double.parseDouble(st.nextToken()); 

                double mfe1 = Nmfe/cg;            
                double mfe4 = mfe/bp;
		double mfe5 = mfe4/cg;
		
                p.println(cg+"	"+bp+"	"+mfe1+"  "+mfe4+"	"+mfe5);
                
                //System.out.println(s);
          }
            

          p.close();
          r1.close();

      }

      catch(Exception e)
      {
         System.out.println(e.toString());
       }
   }

}    







