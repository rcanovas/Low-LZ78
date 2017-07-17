## Low-LZ78

The Low-LZ78 project contains the C++11 and C codes associated with the following reference: 
D. Arroyuelo, R. Canovas, G. Navarro, and R. Raman. "LZ78 Compression in Low Main Memory Space". To appear in Proc. SPIRE'17


This work includes the following data structures:

C++11:
- hlz78: A LZ78 data structure that uses a fixed hash table.
- ghlz78: A LZ78 data structure that uses a growing hash table.
- mhlz78: A LZ78 data structure that uses multi-hash tables.

C (32 bits):
- lz78_min: The LZ78 compact representation of Arroyuelo and Navarro ("Space-efficient construction of Lempel-Ziv compressed text indexes". Information and Computation, 209(7):1070–1102, 2011).
- lz78_uc: The LZ78 uncompressd baseline representation used by Arroyuelo and Navarro.
 

## Compile

To be able to compile the codes: 
- Install sdsl-lite. Follow the installation guide here: (https://github.com/simongog/sdsl-lite)
- Modify the location of the sdsl library in the CMakeLists.txt if necessary.
- Need cmake version 3.5 or higher.
- For the C++11 code, create the build folder and run: 
	- cmake ..
	- make
- For the C codes, go to the lz78_c32 folder and run:
	- make

## Methods

### hlz78/ghlz78/mhlz78

-[CompressText] :

	Use: ./CompressText file_name <opt> 
      		<file_name>: Name of the file to be compressed 
		<opt> : 
		  	-o output_name:  String containing the name of the output file. Default file_name
		  	-s sigma: Size of the alphabet of the input text. Default s = 256
		  	-f factor: Integer value indicating the overload factor used  for the Hash table. Default factor = 5%
		  	-d D_bits: Number of bits used to store each collision value. Default d=0 (meaning that it is compute internally)
			-w Index_type. Default = 0
        		    ---+--------------------
 			     0 | HLZ78 using a map for the displacements 
 			     1 | HLZ78 using a hash table and a sublayer for displacements
 			     2 | MHLZ78 using maps for the displacements
 			     3 | MHLZ78 using a hash table and a sublayer for displacements
 			     4 | GHLZ78 using maps for the displacements
 			     5 | GHLZ78 using a hash table and a sublayer for displacements
        		     

          	output:  <output_name>.index_type

		Example: ./CompressText file -w 1 -s 97 -f 10 -d 3
		output:  file.hlz78_hash
        

-[DecompressHLZ] :

	Use: ./DecompressHLZ index_type <opt>
		<index_type>: Name of the compressed index file
		<opt> (index details are requiered): 
			-o output_name:  String containing the name of the output file. Default <index_name>_decom
			-w Index_type. Default = 0
			 # | Index_type
			 ---+--------------------
			  0 | HLZ78 using a map for the displacements
			  1 | HLZ78 using a hash table and a sublayer for displacements
			  2 | MHLZ78 using maps for the displacements
			  3 | MHLZ78 using a hash table and a sublayer for displacements
			  4 | GHLZ78 using maps for the displacements
			  5 | GHLZ78 using a hash table and a sublayer for displacements
     
		output:  output_name file

		Example: ./DecompressHLZ file.hlz78_hash -w 1 -o file   
		output: file (uncompress version of file.hlz78_hash)
		        

### lz78_min/lz78_uc

-[Clz78_(min/uc)] :
	
	Use: ./Clz78_(min/uc) file_name
		file_name: Name of the file to be compressed 

		output:  file_name.lzt file

		Example: ./Clz78_min file 
		output: file.lzt (Compressed version of file)

-[Dlz78_(min/uc)] :
	
	Use: ./Dlz78_(min/uc) file_name
		file_name: Name of the file to be decompressed (file_name.lzt must exist) 

		output:  file_name file

		Example: ./Dlz78_min file 
		output: file (Decompressed version of file.lzt. Note: the original file is overwrite)


For more information please refer to the paper "LZ78 Compression in Low Main Memory Space". To appear in Proc. SPIRE'17"
	
			
Note: These codes assume that the computer used has enough RAM to read and store the complete input.
