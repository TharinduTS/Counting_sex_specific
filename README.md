# Counting_sex_specific

I started with VCF files and converted into tab delimited files with following commands on computecanada
```bash
# Prepare input file
 module load tabix
 bgzip -c file.vcf > file.vcf.gz
 tabix -p vcf file.vcf.gz
#  Now use vcftools to make a tab delimited file:
 module load StdEnv/2020 vcftools/0.1.16
 zcat file.vcf.gz | vcf-to-tab > out.tab
```
Then wrote following python script to count sex chromosome specific looking sites in windows
#**** PLEASE CHECK WORKING DIRECTORY USED BY SCRIPT BEFORE RUNNING *****************

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 13:57:57 2023

@author: Tharindu
"""

import array

print('count different nucleotides by population')
# set path to the current folder
print('Check the working directory before running')
import os

import inspect
import os


# to get the current working directory
directory = "/Users/Tharindu/Library/CloudStorage/OneDrive-McMasterUniversity/for_lab_and_research/Tharindu_on_Mac/lab/python_projects/real_data/"
os. chdir(directory)
print('working in', directory)
print('\n\n\n\n\n')


inputfile = input('enter the file name with full path to your tab file with genotypes\n')


#sex = input('enter sex of each individual in correct order in capitol letters as M of F (eg. MMFF)\n')
#outputfile1 = input('enter the name you want for your output\n')
#proportion = input('enter the proportion here\n')

#print('Proportion is ', proportion, '\n')

print('printing first 10 lines of the file')


def open_file(inputfile):
    try:
        with open(inputfile, 'r') as file:
            content = file.read()
            # Do something with the file content (e.g., print it)
            return content
            print(content, 10)
    except FileNotFoundError:
        print(f"Error: File '{inputfile}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


file_contents = open_file(inputfile)


#reading first line
tab_header=file_contents.split('\n', 1)[0]

# splitting 
samp_list=tab_header.split('\t')

print('your sample list is')
my_samp_list=samp_list[3:len(samp_list)]
print(my_samp_list)

file = open('my_samp_list.txt','w')
for item in my_samp_list:
	file.write(item+"\n")
file.close()

print('\n\n\n\n\n')

print('Your sample list was created as my_samp_list.txt. Please open with excel and add the corresponding genders as M or F in second column and population names in the third column***AND SAVE IT and enter the file name with full path')
file_with_sex_path = input('enter the file name with full path to your tab file with sexes and pop names\n')


data = [[] for _ in range(len(my_samp_list))]
with open(inputfile) as f:
    for line in f.readlines():
        elements = line.split()
        for x in range(len(my_samp_list)):
            data[x].append(elements[x])
            
samplevise_locs=data[3:len(data)]

print(samplevise_locs)

#saving this dataframe with all the data for later
import pandas as pd

full_data_df=pd.DataFrame(samplevise_locs)
full_data_df=full_data_df.T

print(data[3])

# creating needed dataframes



file_with_sex=pd.read_table(file_with_sex_path,header=None)

file_with_sex.columns=['file','sex','pop']

print('Input was')

print(file_with_sex)


# Split DataFrame by sex
males_data = file_with_sex[file_with_sex['sex'] == 'M'] 

print(males_data)

# get a list of names
names=males_data['pop'].unique().tolist()

print('\n\n')

print('available populations for males are\n')
print(names)
print('\n')

valid_response_list=names
pop_to_test=input('\n Please enter a valid population name\n')

while pop_to_test not in valid_response_list:
    pop_to_test=input('\n Please enter a valid population name\n')


males_data = males_data[males_data['pop'] == pop_to_test]

print('using following males')
print('\n\n')
print(males_data)




# Split DataFrame by sex in females
females_data = file_with_sex[file_with_sex['sex'] == 'F'] 

print(females_data)

# get a list of names
names=females_data['pop'].unique().tolist()

print('\n\n')

print('available populations for females are\n')
print(names)
print('\n')

valid_response_list=names
pop_to_test=input('\n Please enter a valid population name\n')

while pop_to_test not in valid_response_list:
    pop_to_test=input('\n Please enter a valid population name\n')

females_data = females_data[females_data['pop'] == pop_to_test]

print('using following females')
print('\n\n')
print(females_data)



file_with_locs=pd.read_csv(inputfile, sep='\t',header=None)
#switch rows and columns
file_with_locs_transposed=file_with_locs.T


#rename columns
file_with_locs_transposed.columns = ['C'+str(i) for i in range(0,file_with_locs_transposed.shape[1])]

#rename first column to file
file_with_locs_transposed = file_with_locs_transposed.rename(columns={'C0': 'file'})

print(file_with_locs_transposed)

males_data_with_locs=males_data.merge(file_with_locs_transposed, how = 'left', on = ['file'])

print(males_data_with_locs)


females_data_with_locs=females_data.merge(file_with_locs_transposed, how = 'left', on = ['file'])

print(females_data_with_locs)



# Get unique values for columns

#create an empty df first
filtered_df = pd.DataFrame(columns=['loc','Col','genotype','males','females'])

import numpy as np

print('\n\n\n\n\n\n\n\n')

print('Here is the summary of genotypes for males and females\n\n')
print('*********************************************************')
print('\n\n\n\n\n')



# for x in range(1, len(males_data_with_locs.axes[1])-3):
#     progress_perc=100*x/(len(males_data_with_locs.axes[1])-3)
#     print(round(int(progress_perc), 2),'  % completed')
#     current_col='C'+str(x)
    
#     male_genotypes=males_data_with_locs.eval(current_col).value_counts(normalize=True)
    
    
#     male_genotypes_df= male_genotypes.to_frame(name='Males').reset_index()
#     male_genotypes_df.columns=['genotype','males']
#     female_genotypes=females_data_with_locs.eval(current_col).value_counts(normalize=True)
#     female_genotypes_df= female_genotypes.to_frame(name='Females').reset_index()
#     female_genotypes_df.columns=['genotype','females']
#     all_data_summary=male_genotypes_df.merge(female_genotypes_df, how = 'left', on = ['genotype'])
#     all_data_summary.insert(0, 'Col', current_col)
#     all_data_summary.insert(0, 'loc', file_with_locs_transposed.iloc[1,x])
#     all_data_summary=all_data_summary.fillna(0)
#     filtered_df=pd.concat([filtered_df, all_data_summary])
    
# # print(filtered_df)


# for x in range(1, len(males_data_with_locs.axes[1])-3):
#     progress_perc=100*x/(len(males_data_with_locs.axes[1])-3)
#     print(round(int(progress_perc), 2),'  % completed')
#     current_col='C'+str(x)
#     print('counting males')
#     genotype, males = np.unique(males_data_with_locs.eval(current_col), return_counts=1)
#     male_genotypes_df = pd.DataFrame({'genotype': genotype, 'males': list(males)}, columns=['genotype', 'males'])
#     print('counting females')
#     genotype_f, females = np.unique(females_data_with_locs.eval(current_col), return_counts=1)
#     female_genotypes_df = pd.DataFrame({'genotype': genotype, 'females': list(females)}, columns=['genotype', 'females'])
#     all_data_summary=male_genotypes_df.merge(female_genotypes_df, how = 'left', on = ['genotype'])
#     #frequency
#     print('colculating percentage')
#     all_data_summary['males']=all_data_summary['males']/all_data_summary['males'].sum()
#     all_data_summary['females']=all_data_summary['females']/all_data_summary['females'].sum()
#     print('inserting loc and col')
#     all_data_summary.insert(0, 'Col', current_col)
#     all_data_summary.insert(0, 'loc', file_with_locs_transposed.iloc[1,x])
#     filtered_df=pd.concat([filtered_df, all_data_summary])


for x in range(1, len(males_data_with_locs.axes[1])-3):
    progress_perc=100*x/(len(males_data_with_locs.axes[1])-3)
    print(round(int(progress_perc), 2),'  % completed')
    current_col='C'+str(x)
    print('seperating aolumn for males')
    counting_now_m=males_data_with_locs[current_col]
    print(counting_now_m)
    print('counting males')
    male_genotypes=counting_now_m.value_counts(normalize=True)
    print('converting to df')
    male_genotypes_df= male_genotypes.to_frame(name='Males').reset_index()
    male_genotypes_df.columns=['genotype','males']
    print('seperating column for females')
    counting_now_f=females_data_with_locs[current_col]
    print(counting_now_f)
    print('counting females')
    female_genotypes=counting_now_f.value_counts(normalize=True)
    print('converting to df')
    female_genotypes_df= female_genotypes.to_frame(name='Females').reset_index()
    female_genotypes_df.columns=['genotype','females']
    all_data_summary=male_genotypes_df.merge(female_genotypes_df, how = 'outer', on = ['genotype'])
    all_data_summary.insert(0, 'Col', current_col)
    all_data_summary.insert(0, 'loc', file_with_locs_transposed.iloc[1,x])
    all_data_summary=all_data_summary.fillna(0)
    print(all_data_summary)
    filtered_df=pd.concat([filtered_df, all_data_summary])
    

filtered_df=filtered_df.fillna(0)

#going to save this dataframe to be used again as it takes time
saved_big_summary_df = filtered_df

#to loop through this data set
re_filter_dataset = 'y'
    
while (re_filter_dataset == 'y'):
    filtered_df=saved_big_summary_df

    
    
    filtered_df=filtered_df.fillna(0)
    
    print(filtered_df)    
    print('\n\n\n')
    
    print('Do you want to further filter this data?(y/n)\n')
    print('\n\n')
    
    
    valid_response_list='y','n'
    decision = input('enter y or n \n')
    while decision not in valid_response_list:
        decision=input('\n Please enter a valid input\n')
    

    
    if (decision == 'n'):
        print('Thank you. Exiting the script')
        import sys
        sys.exit()
    else:
        print('greater than or less than?')
        
    valid_response_list='g','l'
    decision = input('\n Please enter g or l \n')
    while decision not in valid_response_list:
        decision=input('\n Please enter a valid input\n')
    
    
    print('which sex are you considering?')
    
    valid_response_list='males','females'
    sex_to_check = input('enter males or females \n')
    while sex_to_check not in valid_response_list:
        sex_to_check = input('enter males or females \n')
        
    
    if (decision == 'g'):
        print('enter the value (between 0 and 1 as percentage')
        value_to_check = input('enter value \n')
        while float(value_to_check) > 1:
            value_to_check = input('Please enter a value between 0 and 1 \n')
        value_to_check = float(value_to_check)
        filtered_df=filtered_df.loc[filtered_df[sex_to_check]>=value_to_check]
        print('printing rows with values greater than',value_to_check,'for',sex_to_check)
        print('\n\n')
        print(filtered_df)
    elif (decision == 'l'):
        print('enter the value (between 0 and 1 as percentage')
        value_to_check = input('enter value \n')
        while float(value_to_check) > 1:
            value_to_check = input('Please enter a value between 0 and 1 \n')
        value_to_check = float(value_to_check)
        filtered_df=filtered_df.loc[filtered_df[sex_to_check]<=value_to_check]
        print('printing rows with values less than',value_to_check,'for',sex_to_check)
        print('\n\n')
        print(filtered_df)
    else:
        print('your input was not one of the 2 valid responses')
    
    
    
    
    
    #def countDis(str):
     
        # Stores all distinct characters
        #s = set(str)
     
        # Return the size of the set
       # return len(s)
    
    print('Do you want to filter homozygous or heterozygous sites?\n\n')
    
    valid_response_list='y','n'
    proceed = input('enter y or n \n')
    while proceed not in valid_response_list:
        proceed=input('\n Please enter a valid input\n')
        
    print('\n\n***************************\n')
    print('Preparing data ....')
    print('\n\n***************************\n')
    
    if (proceed=='n'):
        print('Thank you. Exiting the script')
        import sys
        sys.exit()
    else:
        count_list=[]
    
        for i in range(0,len(filtered_df)):
            print(str(round(int(i/len(filtered_df)*100)))+' % completed')
            gen_str=str(filtered_df.iloc[i,2])
            splitted_txt=gen_str.split('/')
            if (splitted_txt[0]==splitted_txt[1]):
                count=2
            else:
                count=3
            # to convert number to list of integers
            next_val = list(map(int, str(count)))
            count_list=count_list+next_val
            
    filtered_df['zygosity'] = count_list
    
    print('homozygous or heterozygous?\n\n')
    
    valid_response_list='hom','het'
    zygo = input('enter hom or het \n')
    while zygo not in valid_response_list:
        zygo=input('\n Please enter a valid input\n')
        
    
    checking='zygosity'
    
    if (zygo=='hom'):
        filtered_df=filtered_df.loc[filtered_df.zygosity==2]
    elif (zygo=='het'): 
        filtered_df=filtered_df.loc[filtered_df.zygosity==3]     
    else:
        print('wrong input detected\n')
    
    print(filtered_df)
    
    
    #saving as tsv file
    
    # exist or not.
    if not os.path.exists("outputs"):
          
        # if the demo_folder directory is not present 
        # then create it.
        os.makedirs("outputs")
    # I have to save this and re import it as it wont work without reading the col name loc 
    
    saving_file_name_for_summary='outputs/'+pop_to_test+'_'+zygo+'_'+sex_to_check+'_'+decision+'_'+str(value_to_check)+'_'+'_filtered_summary.tsv'
    
    filtered_df.to_csv(saving_file_name_for_summary, sep="\t",index=False)
    
    
    # First section ends here*********************
    print('Do you want to do windowed analysis')
    redo_windows='y'
    
    valid_response_list='y','n'
    while redo_windows not in valid_response_list:
        redo_windows=input('\n y or n\n')
        
    keep='y'
    while redo_windows in keep:
        data_to_average=pd.read_table(saving_file_name_for_summary,skiprows=1,header=None)
        data_to_average.columns=('location','Col','genotype','males','females','zygosity')
        #check for window size
        win_size = input('enter_window_size \n')
        win_size=int(win_size)
        step_size= input('enter_step_size \n')
        step_size=int(step_size)
        
        #create an empty dataframe to store windowed data
        windowed_data = pd.DataFrame(columns=['location','number'])
        
        x=1
        jump=win_size
        #count the selected rows between locations
        
        for i in range(1,round(int(max(data_to_average['location'])/win_size))+1):
            sum_of_rows=data_to_average[data_to_average['location'].between(x,jump)].shape[0]
            new_row=x,sum_of_rows
            new_row_df=pd.DataFrame(new_row)
            new_row_df=new_row_df.T
            new_row_df.columns=('location','number')
            windowed_data=pd.concat([windowed_data,new_row_df])
            jump=jump+step_size
            x=x+step_size
            progress_perc=100*i/round(int(max(data_to_average['location'])/win_size))
            print(round(int(progress_perc), 2),'  % completed')
        
        print(windowed_data)
        
        saving_file_name_for_windowed_summary=saving_file_name_for_summary+'_'+str(win_size)+'_'+str(step_size)+'_windowed_summary.tsv'
        
        windowed_data.to_csv(saving_file_name_for_windowed_summary,sep="\t",index=False)
        
        print('\n\nYour output is saved in the directory  outputs')
        
        print('\n\nDo you want to repeat windowed analysis with a different size??')
        redo_windows = input('enter y or n \n')

print('Do you want to do a different filteraing with the same populations before exiting?')
valid_response_list='y','n'
re_filter_dataset = input('enter y or n \n')
while re_filter_dataset not in valid_response_list:
    re_filter_dataset=input('\n Please enter a valid input\n')
```
This creates tsvs with all the filtering methods used in previous script.
First add y w or z according to the filtering steps you used so the next R script can recognize it.
as an example, if your have a file with 'filename' that contains y specific looking sites, the filename should be y_filename

# plotting

Then I used following R script to plot output.

PLEASE CHECK FILE AND POPULATION NAMES BEFORE RUNNING

```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(ggplot2)
library(dplyr)
library(data.table)
library(zoo)

# Delete previous plot before starting if it is in the same folder*****************


#remove scientific notation
options(scipen=999)

#******change pop name here and add y,w,z to the begining of corresponding file names
pop_name<-"Ghana"
#change window size here 
win_size<-50
#change moving win size here 
mov_win_size<-50
#*************

y_sp_sites_name<-paste('y_',pop_name,'_het_females_l_0.0__filtered_summary.tsv_5000_5000_windowed_summary.tsv',sep = '')
w_sp_sites_name<-paste('w_',pop_name,'_hom_males_l_0.0__filtered_summary.tsv_5000_5000_windowed_summary.tsv',sep = '')
z_sp_sites_name<-paste('z_',pop_name,'_hom_females_l_0.0__filtered_summary.tsv_5000_5000_windowed_summary.tsv',sep = '')

y_sp_df<-read.table(y_sp_sites_name,header = TRUE)
w_sp_df<-read.table(w_sp_sites_name,header = TRUE)
z_sp_df<-read.table(z_sp_sites_name,header = TRUE)

file_list<-c(y_sp_sites_name,w_sp_sites_name,z_sp_sites_name)

# combine files adding pop and sex from file name
all_df <- do.call(rbind, lapply(file_list, function(x) cbind(read.table(x, skip = 1), file_name=head(strsplit(x, "_")[[1]],n=1))))


colnames(all_df)<-c('Chr_position','Number_of_sites','Specific_to')

all_df$Number_of_sites<-as.numeric(all_df$Number_of_sites)
#doing another moving average here
#************************
all_df<-subset(all_df, select = c("Chr_position","Number_of_sites","Specific_to"))

#genome coverage as sliding window
all_df_with_average<-setDT(all_df)[, .(
  window.start = rollapply(Chr_position, width=win_size, by=win_size, FUN=min, align="left", partial=TRUE),
  window.end = rollapply(Chr_position, width=win_size, by=win_size, FUN=max, align="left", partial=TRUE),
  average = rollapply(Number_of_sites, width=win_size, by=win_size, FUN=mean, align="left", partial=TRUE)
), .(Specific_to)]

colnames(all_df_with_average)<-c('Specific_to','Chr_position','Number_of_sites','average')
#*************************


max_averages<-all_df_with_average %>% 
  group_by(Specific_to) %>% 
  filter(average==max(average))

colnames(max_averages)<-c('sp','a','b','max')

max_y<-max_averages$max[1]
max_w<-max_averages$max[2]
max_z<-max_averages$max[3]


 #converting to numeric to standardize
 all_df_with_average$Number_of_sites<-as.numeric(all_df_with_average$Number_of_sites)
 all_df_with_average$average<-as.numeric(all_df_with_average$average)
 
 # Adding column based on other column:
all_df_with_average<- all_df_with_average  %>%
   select(Specific_to,Chr_position,Number_of_sites,average) %>%
   mutate(standedized_no = case_when(Specific_to == "y" ~ average/max_y,
                                    Specific_to == "w" ~ average/max_w,
                                    Specific_to == "z" ~ average/max_z))




regular_plot<-ggplot(all_df_with_average,aes(x=Chr_position,y=average,fill=Specific_to))+
  theme_classic()+
  geom_col(position = "dodge")+ 
  scale_fill_manual("legend", values = c("w" = "red", "y" = "blue", "z" = "grey"))+
  ylim(0,320)

standerdize_plot<-ggplot(all_df_with_average,aes(x=Chr_position,y=standedized_no,fill=Specific_to))+
  theme_classic()+
  geom_col(position = "dodge")+ 
  scale_fill_manual("legend", values = c("w" = "red", "y" = "blue", "z" = "grey"))

regular_plot
ggsave(paste(pop_name,"_regular_plot.pdf",sep = ""),width = 15, height = 10)
standerdize_plot
ggsave(paste(pop_name,"_standerdized_plot.pdf",sep = ""),width = 15, height = 10)


```
    
    
    
