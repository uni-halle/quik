package de.uni_halle.barcode_calling.util;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;

public abstract class UtilMethods {
	
	@FunctionalInterface
	public interface TriFunction<A,B,C,R> {
	    R apply(A a, B b, C c);
	}
	
	@FunctionalInterface
	public interface IntBiFunction {
		int apply(int i1, int i2);
	}
	
	public static int min(int... nums) {
		return Arrays.stream(nums).min().orElse(Integer.MAX_VALUE);
	}
	
	public static int max(int... nums) {
		return Arrays.stream(nums).max().orElse(Integer.MAX_VALUE);
	}
	
	public static String getTimestampString() {
		String datePattern = "yyyy-MM-dd_HH-mm-ss";
        DateTimeFormatter df = DateTimeFormatter.ofPattern(datePattern);
            
        //Convert from LocalDateTime to String
        LocalDateTime now = LocalDateTime.now();
        String strNow = df.format(now);
        
        return strNow;
	}
	
	public static int indexOf(int[] array, int val) {
		for (int i = 0; i < array.length; i++)
			if (array[i] == val)
				return i;
		
		return -1;
	}
}
