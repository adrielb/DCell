import javax.comm.*;
import java.io.*;

public class SerialControl implements SerialPortEventListener
{
	private SerialPort serial;
	 OutputStream os;
	 //BufferedInputStream is;
	 InputStream is;
	 char response;
	private StringBuffer availableData = new StringBuffer();
	private byte[] readBuffer = null;
	private String s;
	private int numCommands, numRecieved;
	
	SerialControl(String name, String COM, char resp)
	{
		this( name, COM, resp, 9600);
	}
	
	SerialControl(final String name, String COM, char resp, int baud)
	{
		try
		{
			System.out.println("Initializing: "+name+" on "+COM);
			response=resp;
			CommPortIdentifier port = CommPortIdentifier.getPortIdentifier(COM);
			serial = (SerialPort) port.open(name, 100);
			serial.setSerialPortParams(baud,
			                         SerialPort.DATABITS_8,
			                         SerialPort.STOPBITS_1,
			                         SerialPort.PARITY_NONE);
			serial.setFlowControlMode(serial.FLOWCONTROL_NONE);
			os = serial.getOutputStream();
			//is = new BufferedInputStream(serial.getInputStream());
			is = serial.getInputStream();
			
			serial.addEventListener(this);
			serial.notifyOnDataAvailable(true);
			
			Runtime.getRuntime().addShutdownHook(
				new Thread()
				{
					public void run()
					{
						close();
					}
				});
		} catch (Exception e)
		{
			System.out.println(e);
			System.exit(3);
		}
	}
	
	public synchronized void serialEvent(SerialPortEvent ev)
	{
		try
		{
//System.out.println("Received");
			readBuffer = new byte[is.available()];
			is.read(readBuffer);
			s = new String(readBuffer);
			availableData.append(s);
			
			if(s.indexOf(response)>=0) 
				numRecieved++;
			if( numRecieved==numCommands )
				notify();
//System.out.println("numRecieved: "+numRecieved+"   numC: "+numCommands );
		} catch (Exception e)
		{
			e.printStackTrace(System.out);
		}
	}

	
	public synchronized void write( String command )
	{
		write(command,1);
	}
	
	public synchronized void write( String command, int numCmds )
	{
		try
		{
//System.out.println("Cmd: "+command.trim());
			numCommands = numCmds;
			numRecieved = 0;
			availableData.delete(0,availableData.length());
			os.write( command.getBytes() );
			os.flush();
//System.out.println("Waiting");
			if(response!='!')
				wait();
		} catch (Exception e)
		{
			System.out.println("Error writing: " + command);
			e.printStackTrace(System.out);
		}
	}
	
	
	
	public String getAvailableData()
	{
		return availableData.toString();
	}
	
	public void close()
	{
		if(serial!=null)
		{
			try
			{
				serial.close();
				os.close();
				is.close();
				System.out.println(serial.getName()+" Closed");
			} catch( Exception e )
			{
				e.printStackTrace(System.out);
			}
		}
	}
	
	public static void main(String args[]) throws Exception
	{
		SerialControl sc = new SerialControl("FocusControl","COM5",'\r');
		Thread.sleep(1000);
		long start = System.currentTimeMillis();
		for( int i = 0; i < 1000; i++ )
		{
			sc.write("PZ\r");
//System.out.println("Avail: "+sc.getAvailableData());
		}
		long dur = System.currentTimeMillis()-start;
		System.out.println(dur);
		Thread.sleep(1000);
		sc.close();
	}	
}
