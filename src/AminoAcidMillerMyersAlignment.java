import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static utils.DisplayUtils.printMatrix;

public class AminoAcidMillerMyersAlignment {

    private static final String VALUES_DELIMITER = ",";
    private static final String SEQUENCES_DELIMITER = ":";
    private static final int GAP_PENALTY = -10;
    private final Map<String, Integer> costMatrixMap = new HashMap<>();
    private final Set<String> localOptimalBacktraces = new HashSet<>();
    private int[] leftCostsColumn;
    private int[] rightCostsColumn;
    private int[] leftPrevCostsColumn;
    private int[] rightPrevCostsColumn;



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
                tempSequencesArray = line.split(SEQUENCES_DELIMITER);
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
    }

    public void countAlignmentCost(String[] sequencesArray) {
        String firstSequence = sequencesArray[0];
        String secondSequence = sequencesArray[1];
        System.out.println("Initial sequences:");
        System.out.println(firstSequence);
        System.out.println(secondSequence);
        System.out.println();
        char[] firstSequenceEntries = firstSequence.toCharArray();
        char[] secondSequenceEntries = secondSequence.toCharArray();
        int n = firstSequenceEntries.length;
        int m = secondSequenceEntries.length;
        leftCostsColumn = new int[m + 1];
        rightCostsColumn = new int[m + 1];
        leftPrevCostsColumn = new int[m + 1];
        rightPrevCostsColumn = new int[m + 1];
        int alignmentCost = countLocalAlignmentBestCost(firstSequenceEntries, secondSequenceEntries);
    }

    public int countLocalAlignmentBestCost(char[] dividedSubsequenceEntries, char[] stableSequenceEntries) {
        int n = dividedSubsequenceEntries.length;
        int partitionIndex = n / 2;
        char[] leftSubsequenceEntries = Arrays.copyOfRange(dividedSubsequenceEntries, 0, partitionIndex + 1);
        char[] rightSubsequenceEntries = Arrays.copyOfRange(dividedSubsequenceEntries, partitionIndex, n);
        if (partitionIndex > 0) {
            countLocalAlignmentBestCost(leftSubsequenceEntries, stableSequenceEntries);
            countLocalAlignmentBestCost(rightSubsequenceEntries, stableSequenceEntries);
        }
        if (n == 1) {
            return -1;
        } else {
            leftPrevCostsColumn = getLeftLocalAlignmentColumn(leftSubsequenceEntries, stableSequenceEntries);
            rightPrevCostsColumn = getRightLocalAlignmentColumn(rightSubsequenceEntries, stableSequenceEntries);
        }
    }

    public int[] getLeftLocalAlignmentColumn(char[] firstSequenceEntries, char[] secondSequenceEntries) {
        int i;
        int j;
        int n = firstSequenceEntries.length;
        int m = secondSequenceEntries.length;
        for (j = 0; j < m + 1; j++) {
            leftCostsColumn[j] = j * GAP_PENALTY;
        }
        String replaceCostKey;
        int replaceCostValue;
        for (i = 1; i < n + 1; i++) {
            leftPrevCostsColumn = leftCostsColumn;
            Arrays.fill(leftPrevCostsColumn, 0);
            for (j = 0; j < m + 1; j++) {
                if (j == 0) {
                    leftCostsColumn[j] = i * GAP_PENALTY;
                } else {
                    replaceCostKey = String.valueOf(firstSequenceEntries[i - 1]) + String.valueOf(secondSequenceEntries[j - 1]);
                    replaceCostValue = costMatrixMap.get(replaceCostKey);
                    leftCostsColumn[j] = multiMax(leftPrevCostsColumn[j] + GAP_PENALTY,
                            leftCostsColumn[j - 1] + GAP_PENALTY, leftPrevCostsColumn[j - 1] + replaceCostValue);
                }
            }
        }
        System.out.println(String.format("leftCostsColumn is %s", Arrays.toString(leftCostsColumn)));
    }

    public int[] getRightLocalAlignmentColumn(char[] firstSequenceEntries, char[] secondSequenceEntries) {
        int i;
        int j;
        int n = firstSequenceEntries.length;
        int m = secondSequenceEntries.length;
        int[][] levensteinDistanceMatrix = new int[n + 1][m + 1];
        for (i = n; i >= 0; i--) {
            levensteinDistanceMatrix[i][0] = (n - i) * GAP_PENALTY;
        }
        for (j = m; j >= 0; j--) {
            levensteinDistanceMatrix[0][j] = (n - j) * GAP_PENALTY;
        }
        String replaceCostKey;
        int replaceCostValue;
        for (i = n; i > 0; i--) {
            for (j = m; j > 0; j--) {
                replaceCostKey = String.valueOf(firstSequenceEntries[i - 1]) + String.valueOf(secondSequenceEntries[j - 1]);
                replaceCostValue = costMatrixMap.get(replaceCostKey);
                levensteinDistanceMatrix[i - 1][j - 1] = multiMax(levensteinDistanceMatrix[i][j - 1] + GAP_PENALTY,
                        levensteinDistanceMatrix[i - 1][j] + GAP_PENALTY, levensteinDistanceMatrix[i][j] + replaceCostValue);
            }
        }
        printMatrix(levensteinDistanceMatrix, n + 1, m + 1);
    }

    public int multiMax(int... numbers) {
        int maxNumber = Integer.MIN_VALUE;
        for (int curNumber : numbers) {
            if (curNumber > maxNumber) {
                maxNumber = curNumber;
            }
        }
        return maxNumber;
    }

}
