# include <iostream>
# include <string>
# include <fstream>
# include <vector>

using namespace std;

string path_to_store_uploaded_samples;
string path_single_line_seq_split_prodigal_outputs_tmp;

int number_of_seqs_in_file(string contig_file)
{
	int c = 0;
	string line;
	fstream my_file;
	my_file.open(contig_file, ios::in);
	while (getline (my_file, line)) {
  		
  		if (line[0] == '>'){
  			c++;
  		}
	}
	return c;
}

void create_subfile(string main_fasta_file, string subfile_to_create, int seq_number_to_start_reading, int seq_number_to_finish_reading)
{
	int seq_number = 0;
	string line, seq_description, seq;
	fstream main_contig_file, subfile;
	main_contig_file.open(main_fasta_file, ios::in);
	subfile.open(subfile_to_create, ios::app);
	
	while (getline (main_contig_file, line)) { // first getline() is the seq description starting with '>' and second getline() is the seq itself
  		
  		seq_description = line;
		getline (main_contig_file, line);
		seq = line;
		
		seq_number++;
		if(seq_number >= seq_number_to_start_reading && seq_number <= seq_number_to_finish_reading){
			// create new subfile
			subfile << seq_description + "\n";
			subfile << seq + "\n";
		}
		else{
			if(seq_number > seq_number_to_finish_reading){
				break;
			}
		}
	}
	return;
}

void create_single_line_seq_fasta_file(string main_fasta_file, string new_single_line_fasta_file)
{
	string line, seq_description, seq;
	fstream main_contig_file, new_fasta_file;
	main_contig_file.open(main_fasta_file, ios::in);
	new_fasta_file.open(new_single_line_fasta_file, ios::app);
	
	string prev_seq = "";
	while (getline (main_contig_file, line)) { 
  		
  		if(line[0] == '>'){
  			//previous seq writing
  			
  			if(prev_seq.empty()){
  				// do nothing
  			}
  			else{
  				new_fasta_file << seq_description + "\n";
				new_fasta_file << prev_seq + "\n";
  			}
  			
  			// now this seq
  			seq_description = line;
  			prev_seq = "";
  		}
  		else{
  			prev_seq+=line;
  			
  		}
		
	}
	// handling last sequence in the fasta file
	new_fasta_file << seq_description + "\n";
	new_fasta_file << prev_seq + "\n";
	return;
}


void split_contig_files_for_parallel_prodigal(string contig_file)
{
	int n = 10; // number of chunks to split a file into. Don't choose >26.
	
	int N = number_of_seqs_in_file(contig_file);
	
	if (N >= n){

		int n_seqs_in_each_file = N / n;
		int excess_to_be_put_in_last_file = N % n;
		int total_seqs_in_last_file = n_seqs_in_each_file + excess_to_be_put_in_last_file;
		
		vector<string> id_of_files{"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"};
		int i, first_seq_number_to_start_reading = 1, last_seq_number_to_finish_reading = n_seqs_in_each_file;
		string name_of_new_subfile;
		string path_of_new_subfile;
		for (i=0;i<=n-2;i++)
		{
			name_of_new_subfile = "x" + id_of_files[i] + ".fna"; // subfiles are DNA contigs 
			//path_of_new_subfile = subfile_directory + "/" + name_of_new_subfile;
			create_subfile(contig_file, name_of_new_subfile, first_seq_number_to_start_reading, last_seq_number_to_finish_reading);
			
			first_seq_number_to_start_reading = last_seq_number_to_finish_reading + 1;
			last_seq_number_to_finish_reading = first_seq_number_to_start_reading + n_seqs_in_each_file - 1;
			
		}
		name_of_new_subfile = "x" + id_of_files[i] + ".fna"; // subfiles are DNA contigs
		//path_of_new_subfile = subfile_directory + "/" + name_of_new_subfile;
		
		last_seq_number_to_finish_reading = first_seq_number_to_start_reading + total_seqs_in_last_file - 1;
		create_subfile(contig_file, name_of_new_subfile, first_seq_number_to_start_reading, last_seq_number_to_finish_reading);
	}
	else{ // creating only one subfile

		string name_of_new_subfile = "xa.fna";
		int first_seq_number_to_start_reading = 1;
		int last_seq_number_to_finish_reading = first_seq_number_to_start_reading + N - 1;
		create_subfile(contig_file, name_of_new_subfile, first_seq_number_to_start_reading, last_seq_number_to_finish_reading);
	}
	return;
}

int main(int argc, char *argv[])
{	
	
	string corresponding_contig_file = string(argv[1]);

	string new_contig_file = "new.fna";
	
	create_single_line_seq_fasta_file(corresponding_contig_file, new_contig_file);
	split_contig_files_for_parallel_prodigal(new_contig_file);
	
	return 0;

}
