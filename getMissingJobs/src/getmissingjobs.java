
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class getmissingjobs {
    String parameters[] = new String[2];

    public getmissingjobs (int minDirectory, int maxDirectory, int minPoint, int maxPoint, int maxJobID, String N, String gaps, 
            String G, String deltaDiameter, String initMode, String alpha, String beta, String startArgumentsFileName, String fileRunningJobsName, String filterText, 
            String fileJobIDName, int startPreviousIteration, String[] fileDirectoriesName, String[] filePointsName) {
        
        try { 
            parameters[0]=fileDirectoriesName[1];
            parameters[1]=filePointsName[1];
            
            boolean[][] jobFailedList = new boolean[maxDirectory][maxPoint];
            for (int i=0;i<jobFailedList.length;i++) for (int j=0;j<jobFailedList[0].length;j++) jobFailedList[i][j]=true;
            
            BufferedReader runningJobs = new BufferedReader(new FileReader(fileRunningJobsName));
            String buffer; while ((buffer=runningJobs.readLine())!=null) {
                String result[] = buffer.split(" ");
                for (int i=0;i<result.length;i++) if (result[i].contains(filterText)) {
                    result[i] = result[i].substring(filterText.length());
                    int parametersValues[] = new int[2];
                    for (int j=0;j<parameters.length;j++) {
                        int parameterStartIndex = result[i].indexOf(parameters[j]), parameterEndIndex=parameterStartIndex+1;
                        if (parameterStartIndex>=0) {
                            while (parameterEndIndex<result[i].length() && Character.isDigit(result[i].charAt(parameterEndIndex))) parameterEndIndex++;
                            parametersValues[j]=Integer.parseInt(result[i].substring(parameterStartIndex+1,parameterEndIndex));
                        }
                    }
                    try {jobFailedList[parametersValues[0]-1][parametersValues[1]]=false;} catch (Exception e1) {}  //moga byc parametry WYZSZE [outOfBounds] zadane w innej serii
                    break;
                }
            }
            runningJobs.close();
            
            int licznik=0;
            BufferedReader pressures = new BufferedReader(new FileReader(startArgumentsFileName));
            while ((buffer=pressures.readLine())!=null) licznik++; pressures.close(); 
            String pressureList[] = new String[licznik]; 
            pressures = new BufferedReader(new FileReader(startArgumentsFileName));
            for (int i=0;i<licznik;i++) {
                buffer = pressures.readLine();
                pressureList[i] = buffer.split("\t")[1];
            }
            pressures.close();
            
            BufferedWriter saveDirectories = new BufferedWriter(new FileWriter(fileDirectoriesName[0],false)),
                           savePoints = new BufferedWriter(new FileWriter(filePointsName[0],false)),
                           saveJobIDs = new BufferedWriter(new FileWriter(fileJobIDName,false));
            for (int i=minDirectory-1;i<jobFailedList.length;i++) for (int j=minPoint;j<jobFailedList[0].length;j++) if (jobFailedList[i][j]) {
                saveDirectories.write(String.valueOf(i+1)); saveDirectories.newLine();
                savePoints.write(String.valueOf(j)); savePoints.newLine();

                String actualFolderList[] = new File("3D_N-"+N+"_gaps-"+gaps+"_G-"+G+"_badanie-"+String.valueOf(i+1)+"_D-"+deltaDiameter+"_inM-"+initMode+"_a-"+alpha+"_b-"+beta).list();
                int minJobID = maxJobID;
                if (actualFolderList!=null) for (int k=0;k<actualFolderList.length;k++) 
                    if (actualFolderList[k].contains("arg-"+pressureList[G.equals("1")?j:(pressureList.length-1-j)])) 
                        if (actualFolderList[k].contains("Configurations") && !actualFolderList[k].contains(".err") && !actualFolderList[k].contains("Results.txt") && !actualFolderList[k].contains("transient.txt")) 
                            try{minJobID = Math.min(minJobID,Integer.parseInt(actualFolderList[k].substring(2).split("_")[0]));}catch(Exception e1){}
                saveJobIDs.write(String.valueOf(minJobID+startPreviousIteration)); saveJobIDs.newLine();
            }
            saveDirectories.close(); savePoints.close(); saveJobIDs.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) {       
        new getmissingjobs(Integer.parseInt(args[0]),Integer.parseInt(args[1]),Integer.parseInt(args[2]),Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5],args[6], 
            args[7],args[8],args[9],args[10],args[11],args[12],args[13],args[14],args[15],Integer.parseInt(args[16]),args[17].split("!"),args[18].split("!"));
    }
    
}