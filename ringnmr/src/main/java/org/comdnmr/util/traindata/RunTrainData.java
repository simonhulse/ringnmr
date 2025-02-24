package org.comdnmr.util.traindata;

public class RunTrainData {
    public static void main( String[] args ) throws Exception {
        try {
            DataGenerator generator = new DataGenerator("/home/simonhulse/Documents/projects/ring-guessing/config.json");
        } catch (Exception e) {
            throw e;
        }
    }
}
