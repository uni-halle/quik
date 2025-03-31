package de.uni_halle.barcode_calling.exceptions;

public class MissingParameterException extends RuntimeException {

	private static final long serialVersionUID = -8821098322858074777L;
	
	public MissingParameterException(String message) {
		super(message);
	}

}
