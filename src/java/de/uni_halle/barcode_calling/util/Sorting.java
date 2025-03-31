package de.uni_halle.barcode_calling.util;

import java.util.Arrays;

import de.uni_halle.barcode_calling.util.UtilMethods.IntBiFunction;

public class Sorting {
	
	/**

     * Prevents instantiation.

     */

    private Sorting() {}
    
    /*

     * Tuning parameters.

     */


    /**

     * The maximum number of runs in merge sort.

     */

    private static final int MAX_RUN_COUNT = 67;


    /**

     * The maximum length of run in merge sort.

     */

    private static final int MAX_RUN_LENGTH = 33;


    /**

     * If the length of an array to be sorted is less than this

     * constant, Quicksort is used in preference to merge sort.

     */

    private static final int QUICKSORT_THRESHOLD = 286;


    /**

     * If the length of an array to be sorted is less than this

     * constant, insertion sort is used in preference to Quicksort.

     */

    private static final int INSERTION_SORT_THRESHOLD = 20; // 47;


    /**

     * If the length of a byte array to be sorted is greater than this

     * constant, counting sort is used in preference to insertion sort.

     */

    private static final int COUNTING_SORT_THRESHOLD_FOR_BYTE = 29;


    /**

     * If the length of a short or char array to be sorted is greater

     * than this constant, counting sort is used in preference to Quicksort.

     */

    private static final int COUNTING_SORT_THRESHOLD_FOR_SHORT_OR_CHAR = 3200;
	
    
	/**

     * Sorts the specified range of the array using the given

     * workspace array slice if possible for merging

     *

     * @param a the array to be sorted

     * @param left the index of the first element, inclusive, to be sorted

     * @param right the index of the last element, inclusive, to be sorted

     * @param work a workspace array (slice)

     * @param workBase origin of usable space in work array

     * @param workLen usable size of work array

     */

    public static void sort(int[] a, int left, int right,

                     int[] work, int workBase, int workLen,
                     
                     IntBiFunction comp) {
    	
    	// TODO
    	System.out.println("BIG SORT\n");

        // Use Quicksort on small arrays

        if (right - left < QUICKSORT_THRESHOLD) {

            sort(a, left, right, true, comp);

            return;

        }


        /*

         * Index run[i] is the start of i-th run

         * (ascending or descending sequence).

         */

        int[] run = new int[MAX_RUN_COUNT + 1];

        int count = 0; run[0] = left;


        // Check if the array is nearly sorted

        for (int k = left; k < right; run[count] = k) {

            if (comp.apply(a[k], a[k + 1]) < 0) {  // ascending TODO a[k] < a[k + 1]

                while (++k <= right && comp.apply(a[k - 1], a[k]) <= 0);	// TODO a[k - 1] <= a[k]

            } else if (comp.apply(a[k], a[k + 1]) > 0) {  // descending TODO a[k] > a[k + 1]

                while (++k <= right && comp.apply(a[k - 1], a[k]) >= 0);	// TODO a[k - 1] >= a[k]

                for (int lo = run[count] - 1, hi = k; ++lo < --hi; ) {

                    int t = a[lo]; a[lo] = a[hi]; a[hi] = t;

                }

            } else { // equal

                for (int m = MAX_RUN_LENGTH; ++k <= right && comp.apply(a[k - 1], a[k]) == 0; ) {	// TODO a[k - 1] == a[k]

                    if (--m == 0) {

                        sort(a, left, right, true, comp);

                        return;

                    }

                }

            }


            /*

             * The array is not highly structured,

             * use Quicksort instead of merge sort.

             */

            if (++count == MAX_RUN_COUNT) {

                sort(a, left, right, true, comp);

                return;

            }

        }


        // Check special cases

        // Implementation note: variable "right" is increased by 1.

        if (run[count] == right++) { // The last run contains one element

            run[++count] = right;

        } else if (count == 1) { // The array is already sorted

            return;

        }


        // Determine alternation base for merge

        byte odd = 0;

        for (int n = 1; (n <<= 1) < count; odd ^= 1);


        // Use or create temporary array b for merging

        int[] b;                 // temp array; alternates with a

        int ao, bo;              // array offsets from 'left'

        int blen = right - left; // space needed for b

        if (work == null || workLen < blen || workBase + blen > work.length) {

            work = new int[blen];

            workBase = 0;

        }

        if (odd == 0) {

            System.arraycopy(a, left, work, workBase, blen);

            b = a;

            bo = 0;

            a = work;

            ao = workBase - left;

        } else {

            b = work;

            ao = 0;

            bo = workBase - left;

        }


        // Merging

        for (int last; count > 1; count = last) {

            for (int k = (last = 0) + 2; k <= count; k += 2) {

                int hi = run[k], mi = run[k - 1];

                for (int i = run[k - 2], p = i, q = mi; i < hi; ++i) {

                    if (q >= hi || p < mi && comp.apply(a[p + ao], a[q + ao]) <= 0) {	// TODO a[p + ao] <= a[q + ao]

                        b[i + bo] = a[p++ + ao];

                    } else {

                        b[i + bo] = a[q++ + ao];

                    }

                }

                run[++last] = hi;

            }

            if ((count & 1) != 0) {

                for (int i = right, lo = run[count - 1]; --i >= lo;

                    b[i + bo] = a[i + ao]

                );

                run[++last] = right;

            }

            int[] t = a; a = b; b = t;

            int o = ao; ao = bo; bo = o;

        }

    }


    /**

     * Sorts the specified range of the array by Dual-Pivot Quicksort.

     *

     * @param a the array to be sorted

     * @param left the index of the first element, inclusive, to be sorted

     * @param right the index of the last element, inclusive, to be sorted

     * @param leftmost indicates if this part is the leftmost in the range

     */

    public static void sort(int[] a, int left, int right, boolean leftmost, IntBiFunction comp) {
    	
    	// TODO
    	System.out.println("Sort with");
    	System.out.println("left: " + left);
    	System.out.println("right: " + right);
    	System.out.println("leftmost: " + leftmost);

        int length = right - left + 1;


        // Use insertion sort on tiny arrays

        if (length < INSERTION_SORT_THRESHOLD) {

            if (leftmost) {

                /*

                 * Traditional (without sentinel) insertion sort,

                 * optimized for server VM, is used in case of

                 * the leftmost part.

                 */

                for (int i = left, j = i; i < right; j = ++i) {

                    int ai = a[i + 1];

                    while (comp.apply(ai, a[j]) < 0) {	// TODO ai < a[j]

                        a[j + 1] = a[j];

                        if (j-- == left) {

                            break;

                        }

                    }

                    a[j + 1] = ai;

                }

            } else {

                /*

                 * Skip the longest ascending sequence.

                 */

                do {

                    if (left >= right) {

                        return;

                    }

                } while (comp.apply(a[++left], a[left - 1]) >= 0);	// TODO a[++left] >= a[left - 1]


                /*

                 * Every element from adjoining part plays the role

                 * of sentinel, therefore this allows us to avoid the

                 * left range check on each iteration. Moreover, we use

                 * the more optimized algorithm, so called pair insertion

                 * sort, which is faster (in the context of Quicksort)

                 * than traditional implementation of insertion sort.

                 */

                for (int k = left; ++left <= right; k = ++left) {

                    int a1 = a[k], a2 = a[left];


                    if (comp.apply(a1, a2) < 0) {	// TODO a1 < a2

                        a2 = a1; a1 = a[left];

                    }

                    while (comp.apply(a1, a[--k]) < 0) {	// TODO a1 < a[--k]

                        a[k + 2] = a[k];

                    }

                    a[++k + 1] = a1;


                    while (comp.apply(a2, a[--k]) < 0) {	// TODO a2 < a[--k]

                        a[k + 1] = a[k];
                        
                      //TODO
	                    if (k == 0) {
	                    	System.out.println("k: " + k);
	                    	System.out.println("left: " + left);
	                    	System.out.println("right: " + right);
	                    	System.out.println("leftmost: " + leftmost);
	                    	System.out.println("length: " + length);
	                    }

                    }

                    a[k + 1] = a2;

                }

                int last = a[right];


                while (comp.apply(last, a[--right]) < 0) {	// TODO last < a[--right]

                    a[right + 1] = a[right];

                }

                a[right + 1] = last;

            }

            return;

        }


        // Inexpensive approximation of length / 7

        int seventh = (length >> 3) + (length >> 6) + 1;


        /*

         * Sort five evenly spaced elements around (and including) the

         * center element in the range. These elements will be used for

         * pivot selection as described below. The choice for spacing

         * these elements was empirically determined to work well on

         * a wide variety of inputs.

         */

        int e3 = (left + right) >>> 1; // The midpoint

        int e2 = e3 - seventh;

        int e1 = e2 - seventh;

        int e4 = e3 + seventh;

        int e5 = e4 + seventh;


        // Sort these elements using insertion sort

        if (comp.apply(a[e2], a[e1]) < 0) { int t = a[e2]; a[e2] = a[e1]; a[e1] = t; }	// TODO a[e2] < a[e1]


        if (comp.apply(a[e3], a[e2]) < 0) { int t = a[e3]; a[e3] = a[e2]; a[e2] = t;	// TODO a[e3] < a[e2]

            if (comp.apply(t, a[e1]) < 0) { a[e2] = a[e1]; a[e1] = t; }	// TODO t < a[e1]

        }

        if (comp.apply(a[e4], a[e3]) < 0) { int t = a[e4]; a[e4] = a[e3]; a[e3] = t;	// TODO a[e4] < a[e3]

            if (comp.apply(t, a[e2]) < 0) { a[e3] = a[e2]; a[e2] = t;	// TODO t < a[e2]

            	if (comp.apply(t, a[e1]) < 0) { a[e2] = a[e1]; a[e1] = t; }	// TODO t < a[e1]

            }

        }

        if (comp.apply(a[e5], a[e4]) < 0) { int t = a[e5]; a[e5] = a[e4]; a[e4] = t;	// TODO a[e5] < a[e4]

            if (comp.apply(t, a[e3]) < 0) { a[e4] = a[e3]; a[e3] = t;	// TODO t < a[e3]

            	if (comp.apply(t, a[e2]) < 0) { a[e3] = a[e2]; a[e2] = t;	// TODO t < a[e2]

                	if (comp.apply(t, a[e1]) < 0) { a[e2] = a[e1]; a[e1] = t; }	// TODO t < a[e1]

                }

            }

        }


        // Pointers

        int less  = left;  // The index of the first element of center part

        int great = right; // The index before the first element of right part


        if (comp.apply(a[e1], a[e2]) != 0 && comp.apply(a[e2], a[e3]) != 0 && 
    		comp.apply(a[e3], a[e4]) != 0 && comp.apply(a[e4], a[e5]) != 0) {

            /*

             * Use the second and fourth of the five sorted elements as pivots.

             * These values are inexpensive approximations of the first and

             * second terciles of the array. Note that pivot1 <= pivot2.

             */

            int pivot1 = a[e2];

            int pivot2 = a[e4];


            /*

             * The first and the last elements to be sorted are moved to the

             * locations formerly occupied by the pivots. When partitioning

             * is complete, the pivots are swapped back into their final

             * positions, and excluded from subsequent sorting.

             */

            a[e2] = a[left];

            a[e4] = a[right];


            /*

             * Skip elements, which are less or greater than pivot values.

             */

            while (comp.apply(a[++less], pivot1) < 0);	// TODO a[++less] < pivot1

            while (comp.apply(a[--great], pivot2) < 0);	// TODO a[--great] > pivot2


            /*

             * Partitioning:

             *

             *   left part           center part                   right part

             * +--------------------------------------------------------------+

             * |  < pivot1  |  pivot1 <= && <= pivot2  |    ?    |  > pivot2  |

             * +--------------------------------------------------------------+

             *               ^                          ^       ^

             *               |                          |       |

             *              less                        k     great

             *

             * Invariants:

             *

             *              all in (left, less)   < pivot1

             *    pivot1 <= all in [less, k)     <= pivot2

             *              all in (great, right) > pivot2

             *

             * Pointer k is the first index of ?-part.

             */

            outer:

            for (int k = less - 1; ++k <= great; ) {

                int ak = a[k];

                if (comp.apply(ak, pivot1) < 0) { // Move a[k] to left part TODO ak < pivot1

                    a[k] = a[less];

                    /*

                     * Here and below we use "a[i] = b; i++;" instead

                     * of "a[i++] = b;" due to performance issue.

                     */

                    a[less] = ak;

                    ++less;

                } else if (comp.apply(ak, pivot2) > 0) { // Move a[k] to right part TODO ak > pivot2

                    while (comp.apply(a[great], pivot2) > 0) {	// TODO a[great] > pivot2

                        if (great-- == k) {

                            break outer;

                        }

                    }

                    if (comp.apply(a[great], pivot1) < 0) { // a[great] <= pivot2 TODO a[great] < pivot1

                        a[k] = a[less];

                        a[less] = a[great];

                        ++less;

                    } else { // pivot1 <= a[great] <= pivot2

                        a[k] = a[great];

                    }

                    /*

                     * Here and below we use "a[i] = b; i--;" instead

                     * of "a[i--] = b;" due to performance issue.

                     */

                    a[great] = ak;

                    --great;

                }

            }


            // Swap pivots into their final positions

            a[left]  = a[less  - 1]; a[less  - 1] = pivot1;

            a[right] = a[great + 1]; a[great + 1] = pivot2;


            // Sort left and right parts recursively, excluding known pivots

            sort(a, left, less - 2, leftmost, comp);

            sort(a, great + 2, right, false, comp);


            /*

             * If center part is too large (comprises > 4/7 of the array),

             * swap internal pivot values to ends.

             */

            if (less < e1 && e5 < great) {

                /*

                 * Skip elements, which are equal to pivot values.

                 */

                while (comp.apply(a[less], pivot1) == 0) {	// TODO a[less] == pivot1

                    ++less;

                }


                while (comp.apply(a[great], pivot2) == 0) {	// TODO a[great] == pivot2

                    --great;

                }


                /*

                 * Partitioning:

                 *

                 *   left part         center part                  right part

                 * +----------------------------------------------------------+

                 * | == pivot1 |  pivot1 < && < pivot2  |    ?    | == pivot2 |

                 * +----------------------------------------------------------+

                 *              ^                        ^       ^

                 *              |                        |       |

                 *             less                      k     great

                 *

                 * Invariants:

                 *

                 *              all in (*,  less) == pivot1

                 *     pivot1 < all in [less,  k)  < pivot2

                 *              all in (great, *) == pivot2

                 *

                 * Pointer k is the first index of ?-part.

                 */

                outer:

                for (int k = less - 1; ++k <= great; ) {

                    int ak = a[k];

                    if (comp.apply(ak, pivot1) == 0) { // Move a[k] to left part TODO ak == pivot1

                        a[k] = a[less];

                        a[less] = ak;

                        ++less;

                    } else if (comp.apply(ak, pivot2) == 0) { // Move a[k] to right part TODO ak == pivot2

                        while (comp.apply(a[great], pivot2) == 0) {	// TODO a[great] == pivot2

                            if (great-- == k) {

                                break outer;

                            }

                        }

                        if (comp.apply(a[great], pivot1) == 0) { // a[great] < pivot2 TODO a[great] == pivot1

                            a[k] = a[less];

                            /*

                             * Even though a[great] equals to pivot1, the

                             * assignment a[less] = pivot1 may be incorrect,

                             * if a[great] and pivot1 are floating-point zeros

                             * of different signs. Therefore in float and

                             * double sorting methods we have to use more

                             * accurate assignment a[less] = a[great].

                             */

                            a[less] = pivot1;

                            ++less;

                        } else { // pivot1 < a[great] < pivot2

                            a[k] = a[great];

                        }

                        a[great] = ak;

                        --great;

                    }

                }

            }


            // Sort center part recursively

            sort(a, less, great, false, comp);


        } else { // Partitioning with one pivot

            /*

             * Use the third of the five sorted elements as pivot.

             * This value is inexpensive approximation of the median.

             */

            int pivot = a[e3];


            /*

             * Partitioning degenerates to the traditional 3-way

             * (or "Dutch National Flag") schema:

             *

             *   left part    center part              right part

             * +-------------------------------------------------+

             * |  < pivot  |   == pivot   |     ?    |  > pivot  |

             * +-------------------------------------------------+

             *              ^              ^        ^

             *              |              |        |

             *             less            k      great

             *

             * Invariants:

             *

             *   all in (left, less)   < pivot

             *   all in [less, k)     == pivot

             *   all in (great, right) > pivot

             *

             * Pointer k is the first index of ?-part.

             */

            for (int k = less; k <= great; ++k) {

                if (comp.apply(a[k], pivot) == 0) {	// TODO a[k] == pivot

                    continue;

                }

                int ak = a[k];

                if (comp.apply(ak, pivot) < 0) { // Move a[k] to left part TODO ak < pivot

                    a[k] = a[less];

                    a[less] = ak;

                    ++less;

                } else { // a[k] > pivot - Move a[k] to right part

                    while (comp.apply(a[great], pivot) > 0) {	// TODO a[great] > pivot

                        --great;

                    }

                    if (comp.apply(a[great], pivot) < 0) { // a[great] <= pivot TODO a[great] < pivot

                        a[k] = a[less];

                        a[less] = a[great];

                        ++less;

                    } else { // a[great] == pivot

                        /*

                         * Even though a[great] equals to pivot, the

                         * assignment a[k] = pivot may be incorrect,

                         * if a[great] and pivot are floating-point

                         * zeros of different signs. Therefore in float

                         * and double sorting methods we have to use

                         * more accurate assignment a[k] = a[great].

                         */

                        a[k] = pivot;

                    }

                    a[great] = ak;

                    --great;

                }

            }


            /*

             * Sort left and right parts recursively.

             * All elements from center part are equal

             * and, therefore, already sorted.

             */

            sort(a, left, less - 1, leftmost, comp);

            sort(a, great + 1, right, false, comp);

        }

    }
	
    
    
    
    
    
    public static void quickSort(int arr[], int begin, int end, IntBiFunction comp) {
        if (begin < end) {
            int partitionIndex = partition(arr, begin, end, comp);

            quickSort(arr, begin, partitionIndex-1, comp);
            quickSort(arr, partitionIndex+1, end, comp);
        }
    }

    private static int partition(int arr[], int begin, int end, IntBiFunction comp) {
        int pivot = arr[end];
        int i = (begin-1);

        for (int j = begin; j < end; j++) {
            if (comp.apply(arr[j], pivot) <= 0) {	// TODO arr[j] <= pivot
                i++;

                int swapTemp = arr[i];
                arr[i] = arr[j];
                arr[j] = swapTemp;
            }
        }

        int swapTemp = arr[i+1];
        arr[i+1] = arr[end];
        arr[end] = swapTemp;

        return i+1;
    }
    
    /*
    public static void mergeSort(int[] a, int n, IntBiFunction comp) {
        if (n < 2) {
            return;
        }
        int mid = n / 2;
        int[] l = new int[mid];
        int[] r = new int[n - mid];

        for (int i = 0; i < mid; i++) {
            l[i] = a[i];
        }
        for (int i = mid; i < n; i++) {
            r[i - mid] = a[i];
        }
        mergeSort(l, mid, comp);
        mergeSort(r, n - mid, comp);

        merge(a, l, r, mid, n - mid, comp);
    }
    
    private static void merge(int[] a, int[] l, int[] r, int left, int right, IntBiFunction comp) {
    		 
	    int i = 0, j = 0, k = 0;
	    while (i < left && j < right) {
	        if (comp.apply(l[i], r[j]) <= 0) {	// TODO l[i] <= r[j]
	            a[k++] = l[i++];
	        }
	        else {
	            a[k++] = r[j++];
	        }
	    }
	    while (i < left) {
	        a[k++] = l[i++];
	    }
	    while (j < right) {
	        a[k++] = r[j++];
	    }
	}
	*/
    
    public static void mergeSortIterative(int[] a, int n, IntBiFunction comp) {
    	int[] b = new int[n], temp;
    	int l, r, rest, sizeLeft, sizeRight;
    	
    	for (l = 0; l < n; l += INSERTION_SORT_THRESHOLD) {
    		rest = n - l;
    		r = rest >= INSERTION_SORT_THRESHOLD ? l + INSERTION_SORT_THRESHOLD - 1 : n - 1;
    		
    		insertionSort(a, l, r, comp);
    	}
    	
    	int sizeOld = INSERTION_SORT_THRESHOLD, size = 2 * INSERTION_SORT_THRESHOLD;
    	int counter = 0;
    	
    	while (sizeOld <= n) {
    		for (l = 0; l < n; l += size) {
        		rest = n - l;
        		sizeLeft = rest >= sizeOld ? sizeOld : rest;
        		sizeRight = rest >= size ? sizeOld : rest - sizeOld;
        		
        		merge(b, l, a, l, sizeLeft, a, l + sizeOld, sizeRight, comp);
        	}
        	
        	sizeOld = size;
        	size *= 2;
        	
        	// Swap
        	temp = a;
        	a = b;
        	b = temp;
        	
        	counter++;
    	}
    	
    	// If auxiliary array contains the most recent merge -> copy
    	if (counter / 2 == 1)
    		merge(b, 0, a, 0, n, a, 0, 0, comp);
    }

  //insertion sort to be used once the mergesort partitions become small enough
  public static void insertionSort(int[] a, int left, int right, IntBiFunction comp) {
	  for (int i = left, j = i; i < right; j = ++i) {
          int ai = a[i + 1];
          while (comp.apply(ai, a[j]) < 0) {	// TODO ai < a[j]
              a[j + 1] = a[j];
              if (j-- == left) {
                  break;
              }
          }
          a[j + 1] = ai;
	  }
  }

  //standard merging two sorted half arrays into single sorted array
  static void merge(int[] merged_a, int start_a, int[] half_a1, int start_a1, int size_a1, 
                           int[] half_a2, int start_a2, int size_a2, IntBiFunction comp) {
      int i, j, k;
      int total_s = size_a1 + size_a2;
      for (i = start_a1, j = start_a2, k = start_a; k < (start_a + total_s); k++) {
          // if reached end of first half array, run through the loop 
          // filling in only from the second half array
          if (i >= start_a1 + size_a1) {
              merged_a[k] = half_a2[j++];
          }
          // if reached end of second half array, run through the loop 
          // filling in only from the first half array
          else if (j >= start_a2 + size_a2) {
              merged_a[k] = half_a1[i++];
          }
          // merged array is filled with the smaller element of the two 
          // arrays, in order
          else {
	          merged_a[k] = comp.apply(half_a1[i], half_a2[j]) < 0 ?
	                        half_a1[i++] : half_a2[j++];	// TODO half_a1[i] < half_a2[j]
          }
      }
  }

  //merge sort data during merging without the additional copying back to array
  //all data movement is done during the course of the merges
  static void mergesortNoCopy(int[] a, int[] b, int l, int r, IntBiFunction comp) {
      if (r - l <= INSERTION_SORT_THRESHOLD) {
          insertionSort(a, l, r, comp);
          return;
      }
      int m = (l + r) / 2;
      //switch arrays to msort b thus recursively writing results to b
      mergesortNoCopy(b, a, l, m, comp); //merge sort left
      mergesortNoCopy(b, a, m + 1, r, comp); //merge sort right
      //merge partitions of b into a
      merge(a, l, b, l, m - l + 1, b, m + 1, r - m, comp); //merge
  }

  public static void mergeSort(int[] a, int n, IntBiFunction comp) {
      int[] aux = Arrays.copyOf(a, n);
      mergesortNoCopy(a, aux, 0, n - 1, comp);
  }

}
