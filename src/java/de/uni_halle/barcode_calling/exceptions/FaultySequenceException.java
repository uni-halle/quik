package de.uni_halle.barcode_calling.exceptions;

public class FaultySequenceException extends RuntimeException {

private static final long serialVersionUID = -8821098322858074775L;
	
	public FaultySequenceException(String message) {
		super(message);
	}

}
