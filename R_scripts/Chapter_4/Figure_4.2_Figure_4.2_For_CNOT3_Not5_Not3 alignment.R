###This script was written by Joe and will make plots of individual traces and all traces overlaid.
###Within the parent directory (set below) the traces need to be in a directory named "traces" and there needs to be an empty directory called "plots" to save these too.

#Imports
library(tidyverse)
library(RColorBrewer)
#display.brewer.all()
#library(colorfindr)

date <- "22_08_08"
conditions <- c("_WT", "_Formaldehyde")
dates_and_conditions <- paste0(date, conditions)

#set parent_dir as the folder with the csv file in----
parent_dir <- paste0("\\\\data.beatson.gla.ac.uk/data/JETTLES/sucrose_gradients/", date)

#set the lower and upper limit for the y axis
ylims <- c(0, 1.2)

#make a list of files to plot data from
setwd(file.path(parent_dir, "traces"))
fylenames <- dir(pattern = "*.csv") # creates the list of all the csv files in the directory

#use a for loop to read in the data, reformat it, plot the trace a save the reformatted data to a list
data_list <- list() #creates a list to save the data in
for (fylename in fylenames){
  #read in data
  df <- read_csv(file = fylename, col_names = T, skip = 45) #skip the first 45 rows which contain information about the run settings
  
  #create a continuous scale of cumulative volume to plot against absorbance
  total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
  step <- total_vol / nrow(df)
  df$cum_vol <- seq(step, total_vol, step)
  
  #extract sample name
  sample_name <- str_remove(fylename, "\\.csv") #removes the .csv
  sample_name <- str_replace_all(sample_name, "\\_", " ") #replaces underscores with spaces
  
  df$sample <- rep(sample_name)
  
  #plot data
  trace_plot <- ggplot(data = df, aes(x = cum_vol, y = Absorbance))+
    geom_line(size = 1, color = "red")+
    theme_bw()+
    scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
    ylim(ylims)+
    ylab("absorbance")+
    xlab("volume")+
    ggtitle(sample_name)+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.line = element_line(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_blank())
  
  png(file = file.path(parent_dir, "plots", str_replace(fylename, ".csv", ".png")), height = 400, width = 800)
  print(trace_plot)
  dev.off()
  
  #save data to list
  data_list[[fylename]] <- df
}
all_data <- do.call("rbind", data_list)

trace_plot <- ggplot(data = all_data %>% filter(sample == "Undigested"), aes(x = cum_vol, y = Absorbance))+
  geom_line(size = 1)+
  scale_colour_brewer(palette = "Dark2")+ 
  theme_bw()+
  scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
  ylim(ylims)+
  ylab("absorbance")+
  xlab("volume")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_blank())#,
#legend.position = "bottom")

png(file = file.path(parent_dir, "plots/overlaid_traces.png"), height = 1000, width = 2000, res = 300)
print(trace_plot)
dev.off()


CNOT3 <- unlist(str_split("MADKRKLQGEIDRCLKKVSEGVEQFEDIWQKLHNAANANQKEKYEADLKKEIKKLQRLRDQIKTWVASNEIKDKRQLIDNRKLIETQMERFKVVERETKTKAYSKEGLGLAQKVDPAQKEKEEVGQWLTNTIDTLNMQVDQFESEVESLSVQTRKKKGDKDKQDRIEGLKRHIEKHRYHVRMLETILRMLDNDSILVDAIRKIKDDVEYYVDSSQDPDFEENEFLYDDLDLEDIPQALVATSPPSHSHMEDEIFNQSSSTPTSTTSSSPIPPSPANCTTENSEDDKKRGRSTDSEVSQSPAKNGSKPVHSNQHPQSPAVPPTYPSGPPPAASALSTTPGNNGVPAPAAPPSALGPKASPAPSHNSGTPAPYAQAVAPPAPSGPSTTQPRPPSVQPSGGGGGGSGGGGSSSSSNSSAGGGAGKQNGATSYSSVVADSPAEVALSSSGGNNASSQALGPPSGPHNPPPSTSKEPSAAAPTGAGGVAPGSGNNSGGPSLLVPLPVNPPSSPTPSFSDAKAAGALLNGPPQFSTAPEIKAPEPLSSLKSMAERAAISSGIEDPVPTLHLTERDIILSSTSAPPASAQPPLQLSEVNIPLSLGVCPLGPVPLTKEQLYQQAMEEAAWHHMPHPSDSERIRQYLPRNPCPTPPYHHQMPPPHSDTVEFYQRLSTETLFFIFYYLEGTKAQYLAAKALKKQSWRFHTKYMMWFQRHEEPKTITDEFEQGTYIYFDYEKWGQRKKEGFTFEYRYLEDRDLQ", ""))

Not5 <- unlist(str_split("MSQRKLQQDIDKLLKKVKEGIEDFDDIYEKFQSTDPSNSSHREKLESDLKREIKKLQKHRDQIKTWLSKEDVKDKQSVLMTNRRLIENGMERFKSVEKLMKTKQFSKEALTNPDIIKDPKELKKRDQVLFIHDCLDELQKQLEQYEAQENEEQTERHEFHIANLENILKKLQNNEMDPEPVEEFQDDIKYYVENNDDPDFIEYDTIYEDMGCEIQPSSSNNEAPKEGNNQTSLSSIRSSKKQERSPKKKAPQRDVSISDRATTPIAPGVESASQSISSTPTPVSTDTPLHTVKDDSIKFDNSTLGTPTTHVSMKKKESENDSEQQLNFPPDRTDEIRKTIQHDVETNAAFQNPLFNDELKYWLDSKRYLMQPLQEMSPKMVSQLESSLLNCPDSLDADSPCLYTKPLSLPHPTSIFFPNEPIRFVYPYDVPLNLTNNENDTDNKFGKDSKAKSKKDDDIYSRTSLARIFMKFDLDTLFFIFYHYQGSYEQFLAARELFKNRNWLFNKVDRCWYYKEIEKLPPGMGKSEEESWRYFDYKKSWLARRCGNDFVYNEEDFEKL", ""))

Not3 <- unlist(unlist(str_split("MSNRKLQKEIDIIFKKIQEGLQDFDYHYERYESIQNTEDDSDNQREKEKLANDLKKEIKKLQKFREQIKHWLQNDTVHTLGPVGTSYESKLAENKSLIEDSMETYKLVEKQSKLKTFSNQSIMMTFMDSTNGDDDDDDDDSDDDELFSGSDEDEEDEYSDLPEDAVNAIHFFKDSILQLREQTHKLDHEYEKLAAKKLRKNNLATIEAKKEKIQSTIGNNKFHQKKLKKLLRQLTSGVVTEFNLIWLLKDDLNKYLESNGDYEFTRETELYDDIFNQIAAVEEDYSEIHEDNHDHSHSHLEKTSTANGTHHESATTTATAAKNSSVTTRTSPSPILTKANGHIKELSVESGSPTTSRSKLSKIHTNTSQSSTHSNTSSSTHAQAPRQSPELTSPAIVKTLRPATAPSKPVGTLKWSAAAAAAIPDVQEKPELVETESNGSSVHKTVSKDQYEEKINYQSPAPTQQQQQPQTSVESKTLSTTQSLASSITDNINLDKYKEVIKNSNLTKTELGLFSDMNLARVPPGIQDLVISFASKRNNDNFKVLVDSTEFNQFALPLHKPYLPEIVQPSYYSQYTSTSFKHPTPLVKFQSYWNQVRANFGFQKLMDEIKSLTLQNNADNLPIIAELTFVLFYGFYYGVTPAENLIAESYLFELGWKPYKNQSEDGSNNLSTGSISLSDINGSKSKSNRSSNYYFWFKRVKLISRGEDMGQPSNIEFGDYQVFDLSFWEIFVKYGFKFDYNFCQMEPSKTLF", "")))

CNOT3_250 <- paste0(CNOT3[1:250], collapse = "")
Not5_250 <- paste0(Not5[1:250], collapse = "")
Not3_250 <- paste0(Not3[1:250], collapse = "")
