#include "readgauge.h"

void extrInfo_hmc_float(const char * in1, const char * in2, int len1, int len2, hmc_float * dest){
  char tmp[len2-len1-3];
  strncpy(tmp, &in1[len1+3], len2-len1-3);
  tmp[len2-len1-3] = '\0';
  *dest = atof(tmp);
}

// two strings in xlf-info and inverter-info are complicated because there are several vars saved in them
// this is not a beautiful implementation!!
void extrInfo_beta(char * in1, char * in2, int len1, int len2, hmc_float * dest1, hmc_float * dest2, hmc_float * dest3, hmc_float * dest4){
  char tmp[len2-len1-3]; 
  int cutoff;
  strncpy(tmp, &in1[len1+3], len2-len1-3);
  tmp[len2-len1-3] = '\0';
  //every number is saved with 6 digits after the "."
  //find the "."
  cutoff = strchr(tmp,'.')-tmp+1;
  char beta[cutoff + 6];
  strncpy(beta, tmp, cutoff+6);
  beta[cutoff+6] = '\0';
  *dest1 = atof(beta);
  //cut of the part ", kappa = "
  strcpy(tmp, &tmp[cutoff+6+10]);
  cutoff = strchr(tmp,'.')-tmp+1;
  char kappa[cutoff + 6];
  strncpy(kappa, tmp, cutoff+6);
  kappa[cutoff+6] = '\0';
  *dest2 = atof(kappa);
    //cut of the part ", mu = "
  strcpy(tmp, &tmp[cutoff+6+7]);
  cutoff = strchr(tmp,'.')-tmp+1;
  char mu[cutoff + 6];
  strncpy(mu, tmp, cutoff+6);
  mu[cutoff+6] = '\0';
  *dest3 = atof(mu);
    //cut of the part ", c2_rec = "
  strcpy(tmp, &tmp[cutoff+6+11]);
  cutoff = strchr(tmp,'.')-tmp+1;
  char c2_rec[cutoff + 6];
  strncpy(c2_rec, tmp, cutoff+6);
  c2_rec[cutoff+6] = '\0';
  *dest4 = atof(c2_rec);
}

void extrInfo_kappa(char * in1, char * in2, int len1, int len2, hmc_float * dest1, hmc_float * dest2){
  char tmp[len2-len1-3]; 
  int cutoff;
  strncpy(tmp, &in1[len1+3], len2-len1-3);
  tmp[len2-len1-3] = '\0';
  //every number is saved with 6 digits after the "."
  //find the "."
  cutoff = strchr(tmp,'.')-tmp+1;
  char kappa[cutoff + 6];
  strncpy(kappa, tmp, cutoff+6);
  kappa[cutoff+6] = '\0';
  *dest1 = atof(kappa);
   //cut of the part ", mu = "
  strcpy(tmp, &tmp[cutoff+6+7]);
  cutoff = strchr(tmp,'.')-tmp+1;
  char mu[cutoff + 6];
  strncpy(mu, tmp, cutoff+6);
  mu[cutoff+6] = '\0';
  *dest2 = atof(mu);
}

void extrInfo_int(char * in1, char * in2, int len1, int len2, int * dest){
  char tmp[len2-len1-3];
  strncpy(tmp, &in1[len1+3], len2-len1-3);
  tmp[len2-len1-3] = '\0';
  *dest = (int) atoi(tmp);
}

// the \n at the end is overwritten by \0
void extrInfo_char(char * in1, char * in2, int len1, int len2, char * dest){
  char tmp[len2-len1-4];
  strncpy(tmp, &in1[len1+3], len2-len1-3);
  tmp[len2-len1-4] = '\0';
  strcpy(dest, tmp);
}

void trim2(char * buff){
  int i=0, j=0;  
  int len = (int)strlen(buff);  
  while (i != len) {  
   if (buff[i] != '\n' || buff[i] != ' ')  
     buff[j++] = buff[i];  
   i++;  
  }  
  buff[j]=0;  
}


hmc_error get_XLF_infos(char * filename, hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, 
		      hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date ) {
  FILE * reader;
  reader = fopen(filename, "r");
  if (reader != NULL) {
    //there are " " inf front of most of the labels, the last one was added for a different style
    char * tmparray [] = {"plaquette", " trajectory nr", " beta", "kappa", "mu", "c2_rec", " time", " hmcversion", " mubar", " epsilonbar", " date", " plaquette"};
    char tmp1[512];
    while (!feof(reader)) {
      fgets (tmp1, 512, reader);
      trim2(tmp1);
      if(strncmp(tmparray[0], tmp1, strlen(tmparray[0])) == 0) extrInfo_hmc_float(tmp1,  tmparray[0], strlen(tmparray[0]), strlen(tmp1), plaquettevalue);
      if(strncmp(tmparray[1], tmp1, strlen(tmparray[1])) == 0) extrInfo_int(tmp1,  tmparray[1], strlen(tmparray[1]), strlen(tmp1), trajectorynr);
      if(strncmp(tmparray[2], tmp1, strlen(tmparray[2])) == 0) extrInfo_beta(tmp1,  tmparray[2], strlen(tmparray[2]), strlen(tmp1), beta, kappa, mu, c2_rec);
      if(strncmp(tmparray[6], tmp1, strlen(tmparray[6])) == 0) extrInfo_int(tmp1,  tmparray[6], strlen(tmparray[6]), strlen(tmp1), time);
      if(strncmp(tmparray[7], tmp1, strlen(tmparray[7])) == 0) extrInfo_char(tmp1,  tmparray[7], strlen(tmparray[7]), strlen(tmp1), hmcversion);
      if(strncmp(tmparray[8], tmp1, strlen(tmparray[8])) == 0) extrInfo_hmc_float(tmp1,  tmparray[8], strlen(tmparray[8]), strlen(tmp1), mubar);
      if(strncmp(tmparray[9], tmp1, strlen(tmparray[9])) == 0) extrInfo_hmc_float(tmp1,  tmparray[9], strlen(tmparray[9]), strlen(tmp1), epsilonbar);
      if(strncmp(tmparray[10], tmp1, strlen(tmparray[10])) == 0) extrInfo_char(tmp1,  tmparray[10], strlen(tmparray[10]), strlen(tmp1), date);
      if(strncmp(tmparray[11], tmp1, strlen(tmparray[11])) == 0) extrInfo_hmc_float(tmp1,  tmparray[11], strlen(tmparray[11]), strlen(tmp1), plaquettevalue); 
    }
  }
  else {
      printf("\t\tUnable to open %s\n", filename);
      return HMC_FILEERROR;
  }
  return HMC_SUCCESS;
}

hmc_error get_inverter_infos(char * filename, char * solver, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver, 
		      int * time, char * hmcversion, char * date ) {
  FILE * reader;
  reader = fopen(filename, "r");
  if (reader != NULL) {
    //there are " " inf front of most of the labels
    char * tmparray [] = {"solver", " epssq", " noiter", " kappa", "mu", " time", " hmcversion", " date"};
    char tmp1[512];
    while (!feof(reader)) {
      fgets (tmp1, 512, reader);
      if(strncmp(tmparray[0], tmp1, strlen(tmparray[0])) == 0) extrInfo_char(tmp1,  tmparray[0], strlen(tmparray[0]), strlen(tmp1), solver);
      if(strncmp(tmparray[1], tmp1, strlen(tmparray[1])) == 0) extrInfo_hmc_float(tmp1,  tmparray[1], strlen(tmparray[1]), strlen(tmp1), epssq);
      if(strncmp(tmparray[2], tmp1, strlen(tmparray[2])) == 0) extrInfo_int (tmp1,  tmparray[2], strlen(tmparray[2]), strlen(tmp1), noiter);
      if(strncmp(tmparray[3], tmp1, strlen(tmparray[3])) == 0) extrInfo_kappa(tmp1, tmparray[3], strlen(tmparray[3]), strlen(tmp1), kappa_solver, mu_solver);
      if(strncmp(tmparray[5], tmp1, strlen(tmparray[5])) == 0) extrInfo_int (tmp1,  tmparray[5], strlen(tmparray[5]), strlen(tmp1), time);
      if(strncmp(tmparray[6], tmp1, strlen(tmparray[6])) == 0) extrInfo_char(tmp1,  tmparray[6], strlen(tmparray[6]), strlen(tmp1), hmcversion);
      if(strncmp(tmparray[7], tmp1, strlen(tmparray[7])) == 0) extrInfo_char(tmp1,  tmparray[7], strlen(tmparray[7]), strlen(tmp1), date);
    }
  }
  else {
      printf("\t\tUnable to open %s\n", filename);
      return HMC_FILEERROR;
  }
  return HMC_SUCCESS;
}

// from http://www.codecodex.com/wiki/Remove_blanks_from_a_string#C
void trim(char * buff){
  int i=0, j=0;  
  int len = (int)strlen(buff);  
  while (i != len) {  
   if (buff[i] != ' ')  
     buff[j++] = buff[i];  
   i++;  
  }  
  buff[j]=0;  
}

// http://xmlsoft.org/xmlreader.html
// compile with gcc ReadXML.c $(xml2-config --cflags) -Wall $(xml2-config --libs)

void get_XML_info_simple(xmlTextReaderPtr reader, int numbers[6], char * field) {
  xmlChar *name, *value;
  name = xmlTextReaderName(reader);
  int type = xmlTextReaderNodeType(reader);
  /*unsigned */char * cmpr[] = {"field", "precision", "flavours", "lx", "ly", "lz", "lt"};
  //check if the desired info follows
  //sometimes there are additional " " that have to be removed with trim(string)
  if (strcmp((char*)name, cmpr[0]) ==0 && type == 1) {
    xmlTextReaderRead(reader);
    value = xmlTextReaderValue(reader);
    trim((char*) value);
    strcpy(field, (char*)value);
    xmlFree(value);
  }
  if (strcmp((char*)name, cmpr[1]) ==0 && type == 1) {
    xmlTextReaderRead(reader);
    value = xmlTextReaderValue(reader);
    numbers[0] = atoi((char*)value);
    xmlFree(value);
  }
  if (strcmp((char*)name, cmpr[2]) ==0 && type == 1) {
    xmlTextReaderRead(reader);
    value = xmlTextReaderValue(reader);
    numbers[1] = atoi((char*)value);
    xmlFree(value);
  }
  if (strcmp((char*)name, cmpr[3]) ==0 && type == 1) {
    xmlTextReaderRead(reader);
    value = xmlTextReaderValue(reader);
    numbers[2] = atoi((char*)value);
    xmlFree(value);
  }
  if (strcmp((char*)name, cmpr[4]) ==0 && type == 1) {
    xmlTextReaderRead(reader);
    value = xmlTextReaderValue(reader);
    numbers[3] = atoi((char*)value);
    xmlFree(value);
  }
  if (strcmp((char*)name, cmpr[5]) ==0 && type == 1) {
    xmlTextReaderRead(reader);
    value = xmlTextReaderValue(reader);
    numbers[4] = atoi((char*)value);
    xmlFree(value);
  }
  if (strcmp((char*)name, cmpr[6]) ==0 && type == 1) {
    xmlTextReaderRead(reader);
    value = xmlTextReaderValue(reader);
    numbers[5] = atoi((char*)value);
    xmlFree(value);
  }
  xmlFree(name);
}

hmc_error get_XML_infos(char * filename, int * prec, int * lx, int * ly, int * lz, int *lt, int * flavours, char * field_out ) {
  xmlTextReaderPtr reader;
  int ret;
  int tmpArray[6];
  char field[100];
 
  reader = xmlNewTextReaderFilename(filename);
  if (reader != NULL) {
      ret = xmlTextReaderRead(reader);
      while (ret == 1) {
          get_XML_info_simple(reader, tmpArray, field);
          ret = xmlTextReaderRead(reader);
      }
      xmlFreeTextReader(reader);
      if (ret != 0) {
          printf("%s : failed to parse\n", filename);
      }
  } else {
      printf("Unable to open %s\n", filename);
      return HMC_FILEERROR;
  }
  *prec = tmpArray[0], *flavours = tmpArray[1]; *lx = tmpArray[2], *ly = tmpArray[3], *lz = tmpArray[4], *lt = tmpArray[5], 
  strcpy(field_out, field);
  return HMC_SUCCESS;
}


// read in binary file and save it as readable file
// since tmLQCD always saves data with BigEndian one has to be careful

// get XML Infos: file to be read + parameters
hmc_error read_meta_data(char * file, int * lx, int * ly, int * lz, int * lt, int * prec, char * field_out, int * num_entries, 
		 int * flavours, hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date, 
		 char * solvertype, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver,  int * time_solver, char * hmcversion_solver, char * date_solver, int * fermion){
  FILE *fp;
  int MB_flag, ME_flag, msg, rec, status, first, switcher=0;
  char *lime_type; 
  size_t bytes_pad; 
  n_uint64_t nbytes;

  //possible lime_types
  char * lime_types[] = {
    "propagator-type","xlf-info", "inverter-info", "gauge-scidac-checksum-copy", "etmc-propagator-format", 
    "scidac-binary-data", "scidac-checksum", "ildg-format", "ildg-binary-data" };

  //read lime file
  fp = fopen (file, "r");
  LimeReader *r;
  r = limeCreateReader(fp);
  first = 1; msg = 0;
  //go through the lime-entries
  while( (status = limeReaderNextRecord(r)) != LIME_EOF ){
    if( status != LIME_SUCCESS ) { 
      fprintf(stderr, "\t\tlimeReaderNextRecord returned status = %d\n", status);
      return HMC_STDERR;
    }
    if (MB_flag == 1 || first)
    {
	first = 0;
	rec = 0;
	msg++;
    }      
    rec++;
    //read header data
    nbytes    = limeReaderBytes(r);
    lime_type = limeReaderType(r);
    bytes_pad = limeReaderPadBytes(r);
    MB_flag   = limeReaderMBFlag(r);
    ME_flag   = limeReaderMEFlag(r);
    if(strcmp (lime_types[0],lime_type)==0 ) {
      if (switcher == 0){
	printf("\tfile containts fermion informations\n");    
	switcher ++;
      }
      else{
	printf("\tfile containts informations about more than one fermion: %i\n", switcher);    
	switcher ++;
      }
    }
    //!!read the inverter-infos for FIRST fermion infos only!!
    if(strcmp (lime_types[2],lime_type)==0 && switcher == 1 ) {
      printf("\tfound inverter-infos as lime_type %s...\n", lime_type);     
      *fermion = *fermion +1;
      //!!create tmporary file to read in data, this can be done better
      FILE * tmp;
      char tmp_file_name[L_tmpnam];
      tmpnam(tmp_file_name);
      tmp = fopen(tmp_file_name, "w");
      if(tmp == NULL) {
	printf("\t\terror in creating tmp file\n");
	return HMC_FILEERROR;
      }
      char buffer [nbytes];
      limeReaderReadData ((void*) buffer,(n_uint64_t *) &nbytes, r);
      fwrite(buffer, 1, sizeof(buffer), tmp);
      fclose(tmp); 
      
      int err;
      err = get_inverter_infos(tmp_file_name, solvertype, epssq, noiter, kappa_solver, mu_solver, time_solver, hmcversion_solver, date_solver);
      if (err != 0){
	printf("\t\tfailure reading InverterInfos\n");
	return HMC_STDERR;
      }
      else
	printf("\tsuccesfully read InverterInfos\n");
      
      remove(tmp_file_name);	
    }      
    //!!read XLF info, only FIRST fermion is read!!
    if(strcmp (lime_types[1],lime_type)==0 && switcher < 2 ) {
      printf("\tfound XLF-infos as lime_type %s...\n", lime_type);     
      //!!create tmporary file to read in data, this can be done better
      FILE * tmp;
      char tmp_file_name[L_tmpnam];
      tmpnam(tmp_file_name);
      tmp = fopen(tmp_file_name, "w");
      if(tmp == NULL) {
	printf("\t\terror in creating tmp file\n");
	return HMC_FILEERROR;
      }
      char buffer [nbytes];
      limeReaderReadData ((void*) buffer,(n_uint64_t *) &nbytes, r);
      fwrite(buffer, 1, sizeof(buffer), tmp);
      fclose(tmp); 
      
      int err;
      err = get_XLF_infos(tmp_file_name, plaquettevalue, trajectorynr, beta, kappa, mu, c2_rec, time, hmcversion, mubar, epsilonbar, date);
      if (err != 0){
	printf("\t\tfailure reading XLFInfos\n");
	return HMC_STDERR;
      }
      else
	printf("\tsuccesfully read XLFInfos\n");
      
      remove(tmp_file_name);	
    }
    //!!read ildg format (gauge fields) or etmc-propagator-format (fermions), only FIRST fermion is read!!
    if(  (strcmp (lime_types[4],lime_type)==0 || strcmp (lime_types[7],lime_type)==0)  && switcher < 2 ) {
      printf("\tfound XML-infos as lime_type %s..\n", lime_type);
      //!!create tmporary file to read in data, this can be done better
      FILE * tmp;
      char tmp_file_name[L_tmpnam];
      tmpnam(tmp_file_name);
      tmp = fopen(tmp_file_name, "w");
      if(tmp == NULL) {
	printf("\t\terror in creating tmp file\n");
	return HMC_FILEERROR;
      }
      char buffer [nbytes];
      limeReaderReadData ((void*) buffer,(n_uint64_t *) &nbytes, r);
      fwrite(buffer, 1, sizeof(buffer), tmp);
      fclose(tmp); 
      
      int err;
      err = get_XML_infos(tmp_file_name, prec, lx, ly, lz, lt, flavours, field_out );
      if (err != 0){
	printf("\t\tfailure reading XMLInfos\n");
	return HMC_STDERR;
      }
      else
	printf("\tsuccesfully read XMLInfos\n");
      
      // different sizes for fermions or gauge fields
      if(strcmp(field_out, "diracFermion") == 0){
	//latSize sites, 4 dirac indices, Nc colour indices, 2 complex indices
        *num_entries = (int) (*lx)*(*ly)*(*lz)*(*lt)*NC*NSPIN*2;
      }
      else if(strcmp(field_out, "su3gauge") == 0){
	// latSize sites, 4 links, 2 complex indices -> 9 complex numbers per link
        *num_entries = (int) (*lx)*(*ly)*(*lz)*(*lt)*2*4*9;
      }
      else {
	return HMC_STDERR;
      }
      remove(tmp_file_name);	
    }
  }
  limeDestroyReader(r);
  fclose(fp);
  return HMC_SUCCESS;
} 
    
hmc_error read_binary_data_single(char * file, float * numArray, int num_entries, int filelength ){
  printf("\treading binary file %s..\n",file);
  int length = sizeof(float), i; 
  
  // open file and get length
  FILE * in;
  in  = fopen(file, "r+b");
   // get length of file to check input
  fseek(in, 0L, SEEK_END);
  int filelength_check = ftell(in);
  fseek(in, 0L, SEEK_SET);
  printf("\tlength of file:\t\t%i\n", filelength_check);
  int num_entries_check = filelength_check/length; 
  printf("\tnumber of entries:\t%i\n" ,num_entries_check);
  if (filelength_check != filelength && num_entries_check != num_entries){
    printf("\twrong filelength or number of entries!!\n");
    return HMC_STDERR;
  }
  
  // read in bytes from file
  //char buf[filelength], buf2[filelength];  
  char * buf; char * buf2;
  buf = (char*) malloc(filelength*sizeof(char));
  buf2 = (char*) malloc(filelength*sizeof(char)); 
  int cter = 0;
  while(cter < filelength){
   buf[cter] = fgetc(in);
   cter++;
  }
  fclose(in);
  
  //if endian is little, all floats must be reversed
  if(!ENDIAN){
    printf("\tThe ENDIANNESS of the system is little, bytes must be reversed\n");
    for (i=0; i<filelength; i+=length){
      buf2[i] = buf[i+3];
      buf2[i+1] = buf[i+2];
      buf2[i+2] = buf[i+1];
      buf2[i+3] = buf[i];
    }
  }
  else {
    printf("\tThe ENDIANNESS of the system is big, bytes must not be reversed\n");
    for (i=0; i<filelength; i++){
      buf2[i] = buf[i];
    }
  }
  // convert buf2 to floats
  for(i = 0; i<num_entries; i++){
    numArray[i] = *((float*) &buf2[i*length]);
  }
  return HMC_SUCCESS;
}
    
int read_data_single(char * file, float * num_array_single, int num_entries, char * field_out){
  FILE *fp;
  int MB_flag, ME_flag, msg, rec, status, first, cter=0, err;
  char *lime_type; size_t bytes_pad; n_uint64_t nbytes;

  //possible lime_types
  char * lime_types[] = {
    "propagator-type","xlf-info", " inverter-info", "gauge-scidac-checksum-copy", "etmc-propagator-format", 
    "scidac-binary-data", "scidac-checksum", "ildg-format", "ildg-binary-data"
  };

  //read lime file
  fp = fopen (file, "r");
  LimeReader *r;
  r = limeCreateReader(fp);
  first = 1; msg = 0;
  while( (status = limeReaderNextRecord(r)) != LIME_EOF ){
    if( status != LIME_SUCCESS ) { 
      fprintf(stderr, "limeReaderNextRecord returned status = %d\n", 
	      status);
      return HMC_STDERR;
    }
    if (MB_flag == 1 || first)
    {
	first = 0;
	rec = 0;
	msg++;
    }
    rec++;
    //read header data
    nbytes    = limeReaderBytes(r);
    lime_type = limeReaderType(r);
    bytes_pad = limeReaderPadBytes(r);
    MB_flag   = limeReaderMBFlag(r);
    ME_flag   = limeReaderMEFlag(r);
    //!! read data only for the FIRST entry!!
    if( (strcmp (lime_types[5],lime_type)==0 || strcmp (lime_types[8],lime_type)==0 ) && cter <1) {
      cter ++;
      int filelength = nbytes; 
      //!!create tmporary file to read in data
      //!!this has to be changed
      FILE * tmp;
      char tmp_file_name[L_tmpnam];
      tmpnam(tmp_file_name);
      tmp = fopen(tmp_file_name, "w");
      if(tmp == NULL) {
	printf("\terror in creating tmp file\n");
	return HMC_FILEERROR;
      }
      //this cant be "char buffer [nbytes];" because this can be too big
      char * buffer;
      int buffersize = nbytes*sizeof(char);
      buffer = (char*) malloc(buffersize);
      limeReaderReadData ((void*) buffer,(n_uint64_t *) &nbytes, r);
      fwrite(buffer, sizeof(char), nbytes, tmp);
      fclose(tmp); 
      free(buffer);

      err = read_binary_data_single(tmp_file_name, num_array_single, num_entries, filelength );
      if (err != 0){
	printf("\terror in reading binary data\n");
	return HMC_STDERR;
      }
      remove(tmp_file_name);
      }
  }
  limeDestroyReader(r);
  fclose(fp);
  return HMC_SUCCESS;
}
      
int read_binary_data_double(char * file, double * numArray, int num_entries, int filelength ){
  printf("\treading binary file %s..\n",file);
  int i, length = sizeof(double); 
  
  // open file and get length
  FILE * in;
  in  = fopen(file, "r+b");
   // get length of file to check input
  fseek(in, 0L, SEEK_END);
  int filelength_check = ftell(in);
  fseek(in, 0L, SEEK_SET);
  printf("\tlength of file:\t\t%i\n", filelength_check);
  int num_entries_check = filelength_check/length; 
  printf("\tnumber of entries:\t%i\n" ,num_entries_check);
  if (filelength_check != filelength && num_entries_check != num_entries){
    printf("\twrong filelength or number of entries!!\n");
    return HMC_STDERR;
  }
  
  // read in bytes from file
  //char buf[filelength], buf2[filelength];  
  char * buf; char * buf2;
  buf = (char*) malloc(filelength*sizeof(char));
  buf2 = (char*) malloc(filelength*sizeof(char));
  int cter = 0;
  while(cter < filelength){
   buf[cter] = fgetc(in);
   cter++;
  }
  fclose(in);
  
  //if endian is little, all floats must be reversed
  if(!ENDIAN){
    printf("\tThe ENDIANNESS of the system is little, bytes must be reversed\n");
    for (i=0; i<filelength; i+=length){
      buf2[i] = buf[i+7];
      buf2[i+1] = buf[i+6];
      buf2[i+2] = buf[i+5];
      buf2[i+3] = buf[i+4];
      buf2[i+4] = buf[i+3];
      buf2[i+5] = buf[i+2];
      buf2[i+6] = buf[i+1];
      buf2[i+7] = buf[i];
    }
    
  }
  else {
    printf("\tThe ENDIANNESS of the system is big, bytes must not be reversed\n");
    for (i=0; i<filelength; i++){
      buf2[i] = buf[i];
    }
  }
  
  // convert buf2 to doubles
  for(i = 0; i<num_entries; i++){
    numArray[i] = *((double*) &buf2[i*length]);
  }
  return HMC_SUCCESS;
}
    
int read_data_double(char * file, double * num_array_double, int num_entries, char * field_out){
  FILE *fp;
  int MB_flag, ME_flag, msg, rec, status, first, cter = 0, err;
  char *lime_type; size_t bytes_pad; n_uint64_t nbytes;

  //possible lime_types
  char * lime_types[] = {
    "propagator-type","xlf-info", " inverter-info", "gauge-scidac-checksum-copy", "etmc-propagator-format", 
    "scidac-binary-data", "scidac-checksum", "ildg-format", "ildg-binary-data"
  };

  //read lime file
  fp = fopen (file, "r");
  LimeReader *r;
  r = limeCreateReader(fp);
  first = 1; msg = 0;
  while( (status = limeReaderNextRecord(r)) != LIME_EOF ){
    if( status != LIME_SUCCESS ) { 
      fprintf(stderr, "\tlimeReaderNextRecord returned status = %d\n", 
	      status);
      return HMC_STDERR;
    }
    if (MB_flag == 1 || first)
    {
	first = 0;
	rec = 0;
	msg++;
    }
    rec++;
    //read header data
    nbytes    = limeReaderBytes(r);
    lime_type = limeReaderType(r);
    bytes_pad = limeReaderPadBytes(r);
    MB_flag   = limeReaderMBFlag(r);
    ME_flag   = limeReaderMEFlag(r);
    //!! read data only for the FIRST entry!!
    if( (strcmp (lime_types[5],lime_type)==0 || strcmp (lime_types[8],lime_type)==0 ) && cter <1) {
      cter ++;
      int filelength = nbytes; 
      //!!create tmporary file to read in data
      //!!this has to be changed
      FILE * tmp;
      char tmp_file_name[L_tmpnam];
      tmpnam(tmp_file_name);
      tmp = fopen(tmp_file_name, "w");
      if(tmp == NULL) {
	printf("\terror in creating tmp file\n");
	return HMC_FILEERROR;
      }
      //this cant be "char buffer [nbytes];" because this can be too big
      char * buffer;
      int buffersize = nbytes*sizeof(char);
      buffer = (char*) malloc(buffersize);
      limeReaderReadData ((void*) buffer,(n_uint64_t *) &nbytes, r);
      fwrite(buffer, sizeof(char), nbytes, tmp);
      fclose(tmp);
      free(buffer);
      
      err = read_binary_data_double(tmp_file_name, num_array_double, num_entries, filelength );
      if (err != 0){
	printf("\terror in reading binary data\n");
	return HMC_STDERR;
      }
      remove(tmp_file_name);
    }
  }

  limeDestroyReader(r);
  fclose(fp);
  return HMC_SUCCESS;
}

hmc_error read_tmlqcd_file(char * file, 
			   int * lx, int * ly, int * lz, int * lt, int * prec, char * field_out, int * num_entries, int * flavours, 
			   hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date, 
			   char * solvertype, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver, int * time_solver, char * hmcversion_solver, char * date_solver,
			   hmc_float ** array, int * hmc_prec){
 
  int err, fermion = 0;
  FILE * checker;
  checker = fopen(file, "r");
  if(checker == 0){
    printf("\tcould not open sourcefile!\n\n");
    return HMC_FILEERROR;
  }
  else{
    printf("**********************************************************\n");
    printf("reading tmlqcd-file %s..\n\n", file);
  }
  fclose(checker);

  printf("reading Metadata..\n");
  err = read_meta_data(file, lx, ly, lz, lt, prec, field_out, num_entries, flavours, plaquettevalue, trajectorynr, 
		     beta, kappa, mu, c2_rec, time, hmcversion, mubar, epsilonbar, date, 
		     solvertype, epssq, noiter, kappa_solver, mu_solver, time_solver, hmcversion_solver, date_solver, &fermion);
  if (err != 0){
    printf("\terror in read_meta_data:\t%i\n\n", err);
    return HMC_STDERR;
  }
  else {
    printf("\treading XML-file gave:\n\t\tfield type:\t%s\n\t\tprecision:\t%d\n\t\tlx:\t\t%d\n\t\tly:\t\t%d\n\t\tlz:\t\t%d\n\t\tlt:\t\t%d\n\t\tflavours:\t%d\n", field_out, *prec, *lx, *ly, *lz, *lt, *flavours);
    printf("\treading XLF-data gave:\n\t\tplaquette:\t%f\n\t\ttrajectorynr:\t%i\n\t\tbeta:\t\t%f\n\t\tkappa:\t\t%f\n\t\tmu:\t\t%f\n\t\tc2_rec:\t\t%f\n\t\ttime:\t\t%i\n\t\thmc-version:\t%s\n\t\tmubar:\t\t%f\n\t\tepsilonbar:\t%f\n\t\tdate:\t\t%s\n", *plaquettevalue, *trajectorynr, *beta, *kappa, *mu, *c2_rec, *time, hmcversion, *mubar, *epsilonbar, date);
    if(fermion != 0) 
      printf("\treading inverter-data gave:\n\t\tsolvertype:\t%s\n\t\tepssq:\t\t%.30f\n\t\tnoiter:\t\t%i\n\t\tkappa_solver:\t%f\n\t\tmu_solver:\t%f\n\t\ttime_solver:\t%i\n\t\thmc-ver_solver:\t%s\n\t\tdate_solver:\t%s\n", solvertype, *epssq, *noiter, *kappa_solver, *mu_solver, *time_solver, hmcversion_solver, date_solver);
  }
  if(*hmc_prec != *prec){
    printf("\nthe precision of hmc and sourcefile do not match, will not read data!!!\n\n");
    return -2;
  }
  else{
    printf("reading data..\n");
    //!!note: the read-routines were not changed, the array is just set to the values of the num_array`s
    if (*prec == 32){
      printf("\tfound data in single precision\n");
      float * num_array_single;
      num_array_single = (float*) malloc(*num_entries*sizeof(float));
      err = read_data_single(file, num_array_single, *num_entries, field_out);
      if (err != 0){
        printf("\terror in read_data_single:\t%i\n\n", err);
        return HMC_STDERR;
      }
      else
        printf("\tsuccesfully read in data\n");
  
      *array = (hmc_float*) malloc(*num_entries*sizeof(hmc_float));
      int i;
      for (i = 0; i< *num_entries; i++){
	(*array)[i]=num_array_single[i];
      }
      free(num_array_single);
    }
    else if (*prec == 64){
      printf("\tfound data in double precision\n");
      double * num_array_double;
      num_array_double = (double*) malloc(*num_entries*sizeof(double));
      err = read_data_double(file, num_array_double, *num_entries, field_out);
      if (err != 0){
        printf("\terror in read_data_double:\t%i\n\n", err);
        return HMC_STDERR;
      }
      else
        printf("\tsuccesfully read in data\n");
    
      *array = (hmc_float*) malloc(*num_entries*sizeof(hmc_float));
      int i;
      for (i = 0; i< *num_entries; i++){
	(*array)[i]=num_array_double[i];
      }
     
      free(num_array_double);
    }
    else {
      printf("\tcould not determine precision of data\n\n");
      return HMC_STDERR;
    }
    printf("\nsuccesfully read tmlqcd-file %s..\n", file);
    printf("**********************************************************\n\n");
    return HMC_SUCCESS;
  }
}


hmc_error sourcefileparameters::set_defaults(){
  lx_source = 0; 
  ly_source = 0;
  lz_source = 0;
  lt_source = 0;
  prec_source = 0;
  num_entries_source = 0;
  flavours_source = 0;
  trajectorynr_source = 0;
  time_source = 0;
  time_solver_source = 0;
  noiter_source = 0;
  plaquettevalue_source = 0;
  beta_source = 0;
  kappa_source = 0;
  mu_source = 0;
  c2_rec_source = 0;
  mubar_source = 0;
  epsilonbar_source = 0;
  epssq_source = 0;
  kappa_solver_source = 0;
  mu_solver_source = 0;
  return HMC_SUCCESS;
}


hmc_error sourcefileparameters::readsourcefile(char * file, int precision, hmc_float ** array){

  int lx, ly, lz, lt, prec, num_entries, flavours, trajectorynr, time, time_solver, noiter;
  hmc_float plaquettevalue, beta, kappa, mu, c2_rec, mubar, epsilonbar, epssq, kappa_solver, mu_solver;
  char field_out[100]; char hmcversion[50]; char date[50]; char solvertype[50]; char hmcversion_solver[50]; char date_solver[50];
  int err, prec_tmp;

  //CP
  //this was done because i am lazy
  prec_tmp = precision;
  hmc_float * array_tmp;
  
  err = read_tmlqcd_file( file, 
			  &lx, &ly, &lz, &lt, &prec, field_out, &num_entries, &flavours, 
			  &plaquettevalue, &trajectorynr, &beta, &kappa, &mu, &c2_rec, &time, hmcversion, &mubar, &epsilonbar, date, 
			  solvertype, &epssq, &noiter, &kappa_solver, &mu_solver, &time_solver, hmcversion_solver, date_solver,
			  &array_tmp, &prec_tmp);
  if(err!= 0){
    printf("error in reading file: %s\n\n", file);
    return -1;
  } 

  *array = (hmc_float*) malloc(num_entries*sizeof(hmc_float));
  for (int i = 0; i < num_entries; i++) (*array)[i] = array_tmp[i];

  free(array_tmp);

  //set the source parameters
  val_assign_source(&lx_source, lx);
  val_assign_source(&ly_source, ly);
  val_assign_source(&lz_source, lz);
  val_assign_source(&lt_source, lt);
  val_assign_source(&prec_source, prec);
  val_assign_source(&num_entries_source, num_entries);
  val_assign_source(&flavours_source, flavours);
  val_assign_source(&trajectorynr_source, trajectorynr);
  val_assign_source(&time_source, time);
  val_assign_source(&noiter_source, noiter);
  val_assign_source(&time_solver_source, time_solver);
  val_assign_source(&plaquettevalue_source, plaquettevalue);
  val_assign_source(&beta_source, beta);
  val_assign_source(&kappa_source, kappa);
  val_assign_source(&mu_source, mu);
  val_assign_source(&c2_rec_source, c2_rec);
  val_assign_source(&mubar_source, mubar);
  val_assign_source(&epsilonbar_source, epsilonbar);
  val_assign_source(&epssq_source, epssq);
  val_assign_source(&kappa_solver_source, kappa_solver);
  val_assign_source(&mu_solver_source, mu_solver);
  strcpy(field_source, field_out);
  strcpy(hmcversion_source, hmcversion);
  strcpy(date_source, date);
  strcpy(solvertype_source, solvertype);
  strcpy(hmcversion_solver_source, hmcversion_solver);
  strcpy(date_solver_source, date_solver);

  return HMC_SUCCESS;
}


void sourcefileparameters::val_assign_source(int * out, int in) {
  (*out) = in;
  return;
}

void sourcefileparameters::val_assign_source(hmc_float * out, hmc_float in) {
  (*out) = in;
  return;
}
