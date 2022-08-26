import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import static utils.DisplayUtils.printMatrix;

public class AminoacidSequencesAlignment {

    private static final int GAP_PENALTY = -10;
    private static final String VALUES_DELIMITER = ",";
    private static final String SEQUENCES_DELIMETER = ":";
    private final HashMap<String, Integer> costMatrixMap = new HashMap<>();
    private int[][] levensteinDistanceMatrix;
    private char[] firstSequenceEntries;
    private char[] secondSequenceEntries;
    private AlignmentMethodStrategy alignmentMethodStrategy;

    public void fillCostMatrixMap(String CSVMatrixFilePath) {
        try (FileReader fr = new FileReader(CSVMatrixFilePath);
             BufferedReader bfr = new BufferedReader(fr)) {
            String line;
            String[] headerAminoacidArr;
            String[] tempCostArr;
            headerAminoacidArr = bfr.readLine().trim().split(VALUES_DELIMITER);
            int lineCounter = 0;
            int i;
            while ((line = bfr.readLine()) != null) {
                tempCostArr = line.split(VALUES_DELIMITER);
                System.out.println(Arrays.toString(tempCostArr));
                for (i = 1; i < tempCostArr.length; i++) {
                    costMatrixMap.put(headerAminoacidArr[lineCounter] + headerAminoacidArr[i - 1],
                            Integer.valueOf(tempCostArr[i]));
                }
                lineCounter++;
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException ioException) {
            ioException.printStackTrace();
        }
        System.out.println(costMatrixMap);
    }

    public void countAlignmentCostsAndCheckBacktraces(String sequencesFilePath) {
        List<String[]> sequencesList = new ArrayList<>();
        try (FileReader fr = new FileReader(sequencesFilePath);
             BufferedReader bfr = new BufferedReader(fr)) {
            String line;
            String[] tempSequencesArray;
            int i = 0;
            while ((line = bfr.readLine()) != null) {
                tempSequencesArray = line.split(SEQUENCES_DELIMETER);
                for (i = 0; i < tempSequencesArray.length; i++) {
                    tempSequencesArray[i] = (tempSequencesArray[i]).trim();
                }
                sequencesList.add(tempSequencesArray);
            }
        } catch (IOException ioException) {
            ioException.printStackTrace();
        }

        sequencesList.stream().forEach(sequenceArray -> {
            countAlignmentCost(sequenceArray);
            alignmentMethodStrategy.trackBacktrace(levensteinDistanceMatrix, costMatrixMap, firstSequenceEntries,
                    secondSequenceEntries, GAP_PENALTY);
        });

        //int alignmentEffectivity = 0;

        /*for (String[] curSequenceArray : sequencesList) {
            countAlignmentCost(curSequenceArray);
            alignmentEffectivity += alignmentMethodStrategy.trackBacktrace(levensteinDistanceMatrix, costMatrixMap,
            firstSequenceEntries, secondSequenceEntries, GAP_PENALTY);
        }*/

        //System.out.println(String.format("Alignment effectivity is equal to %d", alignmentEffectivity));
    }

    public void countAlignmentCost(String[] sequencesArray) {
        String firstSequence = sequencesArray[0];
        String secondSequence = sequencesArray[1];
        System.out.println("Initial sequences:");
        System.out.println(firstSequence);
        System.out.println(secondSequence);
        System.out.println();
        firstSequenceEntries = firstSequence.toCharArray();
        secondSequenceEntries = secondSequence.toCharArray();
        int i = 0;
        int j = 0;
        int n = firstSequenceEntries.length;
        int m = secondSequenceEntries.length;
        levensteinDistanceMatrix = new int[n + 1][m + 1];
        for (i = 0; i < n + 1; i++) {
            levensteinDistanceMatrix[i][0] = i * GAP_PENALTY;
        }
        for (j = 0; j < m + 1; j++) {
            levensteinDistanceMatrix[0][j] = j * GAP_PENALTY;
        }
        String replaceCostKey;
        int replaceCostValue;
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                replaceCostKey = String.valueOf(firstSequenceEntries[i]) + String.valueOf(secondSequenceEntries[j]);
                replaceCostValue = costMatrixMap.get(replaceCostKey);
                levensteinDistanceMatrix[i + 1][j + 1] = alignmentMethodStrategy.multiMax(levensteinDistanceMatrix[i + 1][j] + GAP_PENALTY,
                        levensteinDistanceMatrix[i][j + 1] + GAP_PENALTY, levensteinDistanceMatrix[i][j] + replaceCostValue);
            }
        }
        printMatrix(levensteinDistanceMatrix, n + 1, m + 1);
        System.out.println("Alignment cost is " + alignmentMethodStrategy.retrieveAlignmentCost(levensteinDistanceMatrix, n, m));
    }

    public int needlmanWunschMultiMax(int... numbers) {
        double maxNumber = Double.NEGATIVE_INFINITY;
        for (int curNumber : numbers) {
            if (curNumber > maxNumber) {
                maxNumber = curNumber;
            }
        }
        return (int) maxNumber;
    }

    public int smithWatermanMultiMax(int... numbers) {
        int max = needlmanWunschMultiMax(numbers);
        return max > 0 ? max : 0;
    }

    public void setAlignmentMethodStrategy(AlignmentMethodStrategy alignmentMethodStrategy) {
        this.alignmentMethodStrategy = alignmentMethodStrategy;
    }
}
