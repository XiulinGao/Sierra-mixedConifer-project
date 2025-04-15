library(tidyverse)
library(data.table)
#library(ggplot2)

data_main = file.path("~/Sierra-mixedConifer-project/data/FATES_initialization_inventories")
fig_main  = file.path("~/Sierra-mixedConifer-project/results")
files = list.files(data_main,pattern=".csv")
inv_desc = gsub(pattern="_initialization",replacement="",x=files)
inv_desc = gsub(pattern=".csv",replacement="",x=inv_desc)
n_file = length(files)
file_path = file.path(data_main,files)


## PFT info
user_pftinfo  = TRUE

# List of PFTs, in case user_pftinfo is TRUE
n            = 0
pftinfo      = list()
n            = n + 1
pftinfo[[n]] = list( id     = 1
                     , key    = "pft1"
                     , short  = "Pine"
                     , desc   = "Pine tree"
                     , colour = "#74A089"
                     , parse  ="P*i*n*e[E*v*g*r*n]"
)
n            = n + 1
pftinfo[[n]] = list( id     = 2
                     , key    = "pft2"
                     , short  = "Cedar"
                     , desc   = "Incense cedar"
                     , colour = "#9A8822"
                     , parse  ="C*e*d*a*r[E*v*g*r*n]"
)       
n            = n + 1
pftinfo[[n]] = list( id     = 3
                     , key    = "pft3"
                     , short  = "Fir"
                     , desc   = "White fir"
                     , colour = "#F8AFA8"
                     , parse  ="F*i*r[E*v*g*r*n]"
) 
n            = n + 1
pftinfo[[n]] = list( id     = 4
                     , key    = "pft4"
                     , short  = "Oak"
                     , desc   = "Black oak"
                     , colour = "#FDDDA0"
                     , parse  ="F*i*r[E*v*g*r*n]"
) #end of PFT list

pftinfo  = do.call(what=rbind,args=lapply(X=pftinfo,FUN=as_tibble,stringsAsFactors=FALSE))
npfts = nrow(pftinfo)

## plot settings
gg_device  = c("png")     
gg_depth   = 600         
gg_ptsz    = 18         
gg_ptszl   = 26
gg_width   = 17.5       
gg_widthn  = 14.5
gg_height  = 8.5         
gg_units   = "in"  

## dbh info
ndbhs = 13
dbh_lwr = c(0,5,10,15,20,30,40,50,60,70,80,90,100)
dbh_upr = c(dbh_lwr[-1],dbh_lwr[ndbhs]+2*max(diff(dbh_lwr)))
dbh_labs = c( paste0("paste(",dbh_lwr[-ndbhs],"-",dbh_upr[-ndbhs],")")
                        , paste0("paste(",dbh_lwr[ ndbhs],"-infinity)"))
dbh_desc = c( paste0("paste(paste(",dbh_lwr[-ndbhs],"<=D*B*H)<",dbh_upr[-ndbhs],"*c*m)")
                          , paste0("paste( D*B*H >=",dbh_lwr[ndbhs],"*c*m)"))
pft_labs     = pftinfo$short
names(pft_labs)=c("1","2","3","4")
dbh_labs     = c(dbh_labs[5:13], "paste(10-20)")
dbh_labs     = parse(text=dbh_labs)
names(dbh_labs) = c("5","6","7","8","9","10","11","12","13","4")

for(n in sequence(n_file)){
  
  path_now = file_path[n]
  desc_now = inv_desc[n]
  df_now = fread(path_now)
  df_now = df_now                                   %>% 
           mutate(PFT = case_when(PFT == "pine" ~ 1,
                                  PFT == "cedar" ~ 2,
                                  PFT == "fir" ~ 3,
                                  PFT == "oak" ~ 4),
                  dbh = case_when( DBH < 5 ~ 1,
                                   DBH >= 5 & DBH < 10 ~ 2,
                                   DBH >= 10 & DBH < 15 ~ 3,
                                   DBH >= 15 & DBH < 20 ~ 4,
                                   DBH >= 20 & DBH < 30 ~ 5,
                                   DBH >= 30 & DBH < 40 ~ 6,
                                   DBH >= 40 & DBH < 50 ~ 7,
                                   DBH >= 50 & DBH < 60 ~ 8,
                                   DBH >= 60 & DBH < 70 ~ 9,
                                   DBH >= 70 & DBH < 80 ~ 10,
                                   DBH >= 80 & DBH < 90 ~ 11,
                                   DBH >= 90 & DBH < 100 ~ 12,
                                   DBH >= 100 ~ 13 )) 
  
   df_now = df_now                              %>% 
            filter(dbh!= 1 & dbh!=2)                
   df_now = df_now                              %>% 
            mutate(dbh = ifelse(dbh==3,4,dbh),
                   nplant = 1)                  %>%
            group_by(PFT,dbh)                   %>% 
            summarize(nplant = sum(nplant))     %>% 
            ungroup()
   df_now$PFT = factor(df_now$PFT,levels=unique(df_now$PFT))
   df_now$dbh = factor(df_now$dbh,levels=sort(unique(df_now$dbh)))
 

  
  
  szpf_split = ggplot(df_now,aes(dbh,nplant)) + geom_bar(stat = "identity") +
    facet_wrap(.~PFT,scales="free",labeller=labeller(PFT=pft_labs)) + 
    scale_x_discrete(labels = dbh_labs ) +
    labs(title=desc_now,
         y=expression("Number of trees"~"("~plant~ha^-1~")"),
         x="Size class(cm)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  f.output = paste0(desc_now,"-split.png")
  ggsave( filename = f.output
          , plot     = szpf_split
          , device   = "png"
          , path     = fig_main
          , width    = gg_widthn*0.8
          , height   = gg_height*0.8
          , units    = gg_units
          , dpi      = gg_depth
  )#end ggsave
  
  
  
  szpf_stack = ggplot(df_now,aes(x=dbh,y=nplant,fill=PFT)) + geom_bar(stat = "identity") +
    #facet_wrap(.~ens_label,scales="free") + 
    scale_x_discrete(labels = dbh_labs ) +
    scale_fill_manual(values=pftinfo$colour,labels=as_labeller(pft_labs))+
    labs(title=desc_now,
         y=expression("Number of trees"~"("~plant~ha^-1~")"),
         x="Size class(cm)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="none")
  
  f.output = paste0(desc_now,"-stack.png")
  ggsave( filename = f.output
          , plot     = szpf_stack
          , device   = "png"
          , path     = fig_main
          , width    = gg_widthn*0.3
          , height   = gg_height*0.4
          , units    = gg_units
          , dpi      = gg_depth
  )#end ggsave
  
  
  
  
}


