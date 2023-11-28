### Load the function to download traits data 
get_fishbase_data <- function(x="Regalecus-glesne"){
  
  Species_name <- x
  Genus <- strsplit(x,"-")[[1]][1]
  Species <- strsplit(x,"-")[[1]][2]
  cat("Species_name = ",  Species_name, "\n")
  
  url2 <- paste("https://www.fishbase.se/summary/",x,".html",sep="")
  c <- htmlParse(getURLContent(url2, followlocation=T,.opts=RCurl::curlOptions(ssl.verifypeer=F)))
  link_list <- getHTMLLinks(c, externalOnly=T, xpQuery="//a/@href", baseURL = docName(c))
  
  if(length(link_list) == 0){
    stop(paste(x, " is not an accepted name in fishbase, check for spelling mistakes and/or synonyms", sep = ""))
  }
  a1<- getNodeSet(c, "//div ")
  a <- getNodeSet(c, "//span ")
  rm(c)
  
  if (length(a)!=0){
    
    List_env1 <- c("Marine","Freshwater","Brackish","marine","freshwater","brackish")
    List_env2 <- c("bathydemersal", "bathypelagic", "benthopelagic","demersal",
                   "pelagic", "pelagic-neritic", "pelagic-oceanic", "reef-associated")
    clim <- c("Tropical","Temperate","Boreal","Subtropical","Deep-water","Polar") 
    
    env1 <- lapply(1:length(a), function(j){
      unlist(sapply(List_env1,function(w){grep(w,xmlValue(a[[j]]))}))
    })
    
    env2 <- lapply(1:length(a), function(j){
      unlist(sapply(List_env2,function(w){grep(w,xmlValue(a[[j]]))}))
    })
    
    env3 <- lapply(1:length(a), function(j){
      unlist(sapply(c(List_env2,List_env1),function(w){grep(w,xmlValue(a[[j]]))}))
    })
    
    resOO <- c(which(lapply(env1,length)!=0), which(lapply(env2,length)!=0))
    if(length(resOO)==2){ cond <- min(resOO)
    } else {
      cond <- as.numeric(names(which(table(resOO)==2)))
      if (length(cond)>1){
        r <- which(sapply(1:length(cond), function (w) grep("depth range",xmlValue(a[[cond[w]]])))==1)
        
        if (length(r)==0) {
          xx <- which.max(lapply(env3,length))
          r <- which(cond==xx)
        }
        
        cond <- cond[r]    
      } 
    }
    
    d <- xmlValue(a[[cond]])
    env1 <-  unlist(sapply(List_env1,function(x)regmatches(d,gregexpr(pattern=x,d))))
    
    if(length(env1)==1){env1 <- env1[[1]]}
    if(length(env1)==2){env1 <- paste(env1[1],env1[2],sep="_")}
    if(length(env1)==3){env1 <- paste(env1[1],env1[2],env1[3],sep="_")}
    
    
    env <-  unlist(sapply(List_env2,function(x)regmatches(d,gregexpr(pattern=x,d))))
    if (length(env)==1){ env2 <- env[[1]]}
    if(length(env)>1){
      if (nchar(env[[1]])>nchar(env[[2]])){
        env2 <- env[[1]]
      } else {
        env2 <- env[[2]]
      }
    }
    
    climate <-  unlist(sapply(clim,function(x)regmatches(d,gregexpr(pattern=x,d))))
    if(length(climate)==0){
      climate <- NA
    } else { 
      climate <- climate[[1]]
    }
    
    if(length(env)==0) {env2 <- NA}
    if(length(env1)==0) {env1 <- NA}
    
    temp <- regmatches(d,gregexpr(pattern= "[0-9]+°C", d))
    if(length(temp[[1]])==0){
      temp_max <- temp_min <-NA
    } else{
      temp_min <-as.numeric(sub(pattern="°C",replacement="",temp[[1]][1]))
      temp_max <-as.numeric(sub(pattern="°C",replacement="",temp[[1]][2]))
    } # end of else
    
    depth_max <- as.numeric(gsub( pattern="[[:space:]]m",replacement="",unlist(regmatches(d,gregexpr(pattern= "[0-9]+ m", d)))[1]))  	
    depth_min <- as.numeric(gsub(pattern="[[:space:]]-",replacement="",unlist(regmatches(d,gregexpr(pattern= "[0-9]+ -", d)))[1]))
    
    depth_max_us <- as.numeric(gsub( pattern="[[:space:]]m",replacement="",unlist(regmatches(d,gregexpr(pattern= "[0-9]+ m", d)))[2]))  	
    depth_min_us <-as.numeric(gsub(pattern="[[:space:]]-",replacement="",unlist(regmatches(d,gregexpr(pattern= "[0-9]+ -", d)))[2]))
    
    b <- xmlValue(a[[which.max(sapply(lapply(a,xmlValue), function(x){regexec(pattern="Max length", x)[[1]][1]}))]])
    rm(a)
    
    Max_length <- substr(b, regexec(pattern="Max length", b)[[1]][1]+13, regexec(pattern="Max length", b)[[1]][1]+22)
    Max_length <-  strsplit(Max_length,split=" ")[[1]][1]
    
    if(sum(grep(pattern=",",Max_length))==1){Max_length = sub (pattern = ",",replacement = "",Max_length) }
    
    if(length(grep("common length",b))==0){
      Common_length = NA
    } else {
      Common_length = substr(b, regexec(pattern="common length", b)[[1]][1]+16, regexec(pattern="common length", b)[[1]][1]+23)
      Common_length = gsub("[[:alpha:]]", "", gsub("[[:space:]]*", "", Common_length, perl=TRUE), perl=TRUE)
      if(sum(grep(pattern=",",Common_length))==1){Common_length = sub(pattern = ",",replacement = "",Common_length) }
      if(Common_length=="y:" | Common_length==":&") {Common_length <- NA}
    }
    
    
    
    ##### Resilience 
    w <- which(sapply(lapply(a1,xmlValue), function(x){regexec(pattern="Resilience", x)[[1]][1]})>0)
    if(length(w)==0){
      Resilience <- NA
    } else {
      d1       <- xmlValue(a1[[w[length(w)]]])
      Res <- regmatches(d1,gregexpr(pattern= "[[:alpha:]]+,", d1))
      Resilience <- sub(pattern=",",replacement="",Res[[1]][1])
    } # end of ifesle
    
    # Vulnérabilité
    w_Vul <- which(sapply(lapply(a1,xmlValue), function(x){regexec(pattern="Vulnerability", x)[[1]][1]})>0)
    if(length(w_Vul)==0){
      Vul <- NA
    } else {
      d1_Vul <- xmlValue(a1[[w_Vul[length(w_Vul)]]])
      Vul <- strsplit(d1_Vul,split=":" )[[1]][2]
      Vul <- regmatches(Vul, gregexpr(pattern= "[[:digit:]]",Vul ))
      if(length(Vul)>5){
        Vul=100
      } else{
        Vul=as.numeric(paste(Vul[[1]][1],Vul[[1]][2],sep=""))
      } # end of ifesle
    } # end of ifesle
    
    ##### Price category
    w_Price <-which(sapply(lapply(a1,xmlValue), function(x){regexec(pattern="Price category", x)[[1]][1]})>0)
    if(length(w_Price)==0){
      Price=NA
    } else {
      d1_Price <- xmlValue(a1[[w_Price[length(w_Price)]]])
      int <- strsplit(d1_Price,":")[[1]][2]
      Price <- unlist(regmatches(int, gregexpr(pattern= "[[:alpha:]]+.",int )))
      if(length(Price)>1){
        Price[1] <- sub(" ","",Price[1])
        Price[2] <- sub(pattern="[[:punct:]]",replacement="",Price[2])
        Price <- paste(Price[1],Price[2],sep="_")
      } else {
        Price <- sub(" ","",Price)
        Price <- sub(pattern="[[:punct:]]",replacement="",Price)
      } # end of ifelse	
    } # end of ifelse
    
    #### Repro
    Repro_link <- link_list[grep("FishReproSummary",link_list)]
    Repro_link <-Repro_link[which.max(nchar(unique(Repro_link)))]
    Repro_link <- gsub("..","",Repro_link,fixed=T)
    
    url_end <- paste ("https://www.fishbase.se/",Repro_link,sep="")
    tt <- getURLContent(url_end, followlocation=TRUE, .encoding="CE_UTF8",.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
    res.t <- readHTMLTable(tt,header=TRUE,colClasses=NULL,skip.rows=integer(),
                           stringsAsFactors=FALSE,trim=TRUE,elFun=xmlValue,
                           as.data.frame=TRUE,which=integer())$dataTable
    
    if(is.null(res.t)){
      Repro.Mode = NA
      Repro.Fertil = NA
      Repro.MatingT = NA
      Repro.SpawnFreq = NA
      Repro.ParentCare = NA
    } else {
      Repro.Mode = res.t[which(res.t[,1] == "Mode"),2] # Mode
      Repro.Fertil = res.t[which(res.t[,1] == "Fertilization"),2] # Fertilization
      Repro.MatingT = res.t[which(res.t[,1] == "Mating type"),2] # Mating type
      Repro.SpawnFreq = res.t[which(res.t[,1] == "Spawning frequency"),2] # Spawning frequency
      Repro.ParentCare = res.t[which(res.t[,1] == "Parental Care"),2] # Parental Care
      
      if(Repro.Mode == "") {Repro.Mode = NA}
      if(Repro.Fertil == "") {Repro.Fertil = NA}
      if(Repro.MatingT == "") {Repro.MatingT = NA}
      if(Repro.SpawnFreq == "") {Repro.SpawnFreq = NA}
      if(Repro.ParentCare == "") {Repro.ParentCare = NA}
    }
    
    #### Eggs
    Eggs_link <- link_list[grep("FishEggInfoSummary",link_list)]
    Eggs_link <-Eggs_link[which.max(nchar(unique(Eggs_link)))]
    Eggs_link <- gsub("..","",Eggs_link,fixed=T)
    
    url_egg <- paste ("https://www.fishbase.se/",Eggs_link,sep="")
    tt <- getURLContent(url_egg, followlocation=TRUE, .encoding="CE_UTF8",.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
    res.t <- tryCatch(readHTMLTable(tt,header=TRUE,colClasses=NULL,skip.rows=integer(),
                                    stringsAsFactors=FALSE,trim=TRUE,elFun=xmlValue,
                                    as.data.frame=TRUE,which=integer())[[1]], error=function(e) paste("Error: Eggs is absent for ", x, sep=""))
    
    if(is(res.t, "character")){
      Egg.Shape = NA
      Egg.Attrib = NA
      Egg.Color = NA
    } else {
      Egg.Shape = res.t[which(res.t[,1] == "Shape of Egg"),2] # Shape of Egg
      Egg.Attrib = res.t[which(res.t[,1] == "Attributes"),2] # Attributes
      Egg.Color = res.t[which(res.t[,1] == "Color of Eggs"),2] # Color of Eggs
      if(Egg.Shape == "") {Egg.Shape = NA}
      if(Egg.Attrib == "") {Egg.Attrib = NA}
      if(Egg.Color == "") {Egg.Color = NA}
    }
    
    #### Larvae
    Larvae_link <- link_list[grep("LarvaeInfoList",link_list)]
    Larvae_link <- Larvae_link[which.max(nchar(unique(Larvae_link)))]
    Larvae_link <- gsub("..","",Larvae_link,fixed=T)
    
    # url_Larvae <- paste ("http://www.fishbase.org/",Larvae_link,sep="")
    url_Larvae <- paste ("https://www.fishbase.se/",Larvae_link,sep="")
    tt <- getURLContent(url_Larvae, followlocation=TRUE, .encoding="CE_UTF8",.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
    res.t <- readHTMLTable(tt,header=TRUE,colClasses=NULL,skip.rows=integer(),
                           stringsAsFactors=FALSE,trim=TRUE,elFun=xmlValue,
                           as.data.frame=TRUE,which=integer())
    
    if(length(res.t) == 0){
      Larvae.LenghtAtBirth = rep(NA, 3)
      names(Larvae.LenghtAtBirth) = c("Larvae.LengthAtBirth.mm.Max", "Larvae.LengthAtBirth.mm.Min", 		"Larvae.LengthAtBirth.mm.Mode")
      Larvae.PlaceDevelop = NA
      names(Larvae.PlaceDevelop) = c("Larvae.PlaceOfDevelopment")
      Larvae.YolkSac = NA
      names(Larvae.YolkSac) = c("Larvae.YolkSac")
    } else {
      # Larvae.LenghtAtBirth
      Larvae.LenghtAtBirth = unlist(lapply(seq(1, length(res.t)), function(x){
        if(length(grep("Length at birth", res.t[[x]])) > 0){
          res.t[[x]][grep("Length at birth", res.t[[x]][,1]),]
        }
      }))[2:4]
      if(is.null(Larvae.LenghtAtBirth)){
        Larvae.LenghtAtBirth = rep(NA, 3)
        names(Larvae.LenghtAtBirth) = c("Larvae.LengthAtBirth.mm.Max", "Larvae.LengthAtBirth.mm.Min", 		"Larvae.LengthAtBirth.mm.Mode")
      } else {
        options(warn=-1)
        Larvae.LenghtAtBirth = as.numeric(as.character(Larvae.LenghtAtBirth))
        options(warn=0)
        names(Larvae.LenghtAtBirth) = c("Larvae.LengthAtBirth.mm.Max", "Larvae.LengthAtBirth.mm.Min", 		"Larvae.LengthAtBirth.mm.Mode")
      }
      
      # Place of development
      Larvae.PlaceDevelop =  unlist(lapply(seq(1, length(res.t)), function(x){
        if(length(grep("Place of development", res.t[[x]])) > 0){
          res.t[[x]][grep("Place of development", res.t[[x]][,1]),]
        }
      }))
      if(is.null(Larvae.PlaceDevelop)){
        Larvae.PlaceDevelop = NA
        names(Larvae.PlaceDevelop) = c("Larvae.PlaceOfDevelopment")
      } else {
        names(Larvae.PlaceDevelop) = c("Larvae.PlaceOfDevelopment")
      }
      
      # Yolk-sac
      Larvae.YolkSac =  unlist(lapply(seq(1, length(res.t)), function(x){
        if(length(grep("Yolk-sac", res.t[[x]])) > 0){
          res.t[[x]][grep("Yolk-sac", res.t[[x]][,1]),2]
        }
      }))[2]
      if(is.null(Larvae.YolkSac)){
        Larvae.YolkSac = NA
        names(Larvae.YolkSac) = c("Larvae.YolkSac")
      } else {
        names(Larvae.YolkSac) = c("Larvae.YolkSac")
      }
    } # End if(length(res.t) == 0){
    
    
    ##### IUCN Status
    w_IUCN  <-which(sapply(lapply(a1,xmlValue), function(x){regexec(pattern="IUCN", x)[[1]][1]})>0)
    if(length(w_IUCN)==0){ 
      IUCN_status=NA
    } else {
      d1_IUCN  <- xmlValue(a1[[w_IUCN[length(w_IUCN)]]])
      IUCN <- unlist(regmatches(d1_IUCN,gregexpr(pattern= "[[:alpha:]]+)", d1_IUCN)))
      IUCN_status <- sub(pattern="[[:punct:]]",replacement="",IUCN[1] ) 
    } # end of ifelse 
    
    ##### Trophic level
    w_TL <- which(sapply(lapply(a1,xmlValue), function(x){regexec(pattern="Trophic level", x)[[1]][1]})>0)
    if(length(w_TL)==0){
      Trophic_Level=NA; TL_var=NA
    } else {
      d1_TL <- xmlValue(a1[[w_TL[length(w_TL)]]])
      TL <- unlist(regmatches(d1_TL,gregexpr(pattern= "[0-9].[0-9]", d1_TL)))
      TL_var <- TL[3] 	
      Trophic_Level <- TL[2]	
    }# end of ifelse
    
    rm(a1) 
    
    ##### Length-weight   
    if(length(grep (pattern = "LW",link_list))==0){
      
      Length_weight = data.frame(a=NA,b=NA,Length_type=NA,n=NA)
      
    } else {
      link_url <-  substr(link_list [[grep (pattern = "LW",link_list)]],2,nchar(link_list [[grep (pattern = "LW",link_list)]]))  
      url_end <- paste ("https://www.fishbase.se/", link_url,sep="")
      
      tt <- getURLContent(url_end, followlocation=TRUE, .encoding="CE_UTF8",.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
      res <- readHTMLTable(tt,header=TRUE,colClasses=NULL,skip.rows=integer(),
                           stringsAsFactors=FALSE,trim=TRUE,elFun=xmlValue,
                           as.data.frame=TRUE,which=integer())[[3]]
      if(is.null(res)==T){
        Length_weight = data.frame(a=NA,b=NA,Length_type=NA,n=NA)
      } else{
        selected_length <- names(which.max(table(res[,"Lengthtype"])))
        
        res <- res[which(res$Lengthtype==selected_length),]
        res$n <- as.character(gsub("\u00A0","",res$n))
        
        if(sum(res$n=="")>=1) {res <- res[-which(res$n==""),]}
        
        n <- sum(as.numeric(as.character(res$n)))
        a <- mean(as.numeric(as.character(res[,"a"])),na.rm=TRUE)
        b <- mean(as.numeric(as.character(res[,"b"])),na.rm=TRUE)
        
        Length_weight <- data.frame(a=a,b=b,Length_type=selected_length,n=n)
        
      }
    }  	
    
    res <- data.frame(c(Max_length= as.numeric(Max_length),Common_length=Common_length,
                        Depth_min=depth_min,Depth_max=depth_max,Depth_min_us=depth_min_us,
                        Depth_max_us=depth_max_us,Env_1=env1,Env_2=env2,Climate=climate,
                        Price=Price,Length_weight,Resilience=Resilience,Vul=Vul,
                        Trophic_Level=Trophic_Level,IUCN_status=IUCN_status,
                        Temp_min=temp_min,Temp_max=temp_max, Repro.Mode=Repro.Mode, Repro.Fertil=Repro.Fertil,
                        Repro.MatingT=Repro.MatingT, Repro.SpawnFreq=Repro.SpawnFreq, 
                        Repro.ParentCare=Repro.ParentCare, Egg.Shape=Egg.Shape, Egg.Attrib=Egg.Attrib,
                        Egg.Color=Egg.Color, Larvae.LenghtAtBirth, Larvae.PlaceDevelop, Larvae.YolkSac))
    
    rownames(res)=x
    
  } else {
    
    Length_weight <- data.frame(a="A_verifier",b="A_verifier",Length_type="A_verifier",n="A_verifier")
    Max_length=Common_length=Resilience=depth_max=depth_min=env1=env2=depth_min_us=Climate=climate=depth_max_us=Price=Trophic_Level=IUCN_status=IUCN_status=temp_min=temp_max=Vul= "A_verifier"
    
    res <- data.frame(c(Max_length= Max_length,Common_length=Common_length,Depth_min=depth_min,Depth_max=depth_max,Depth_min_us=depth_min_us,Depth_max_us=depth_max_us,Env_1=env1,Env_2=env2,Climate=climate,Price=Price,Length_weight,
                        Resilience=Resilience,Vul=Vul,Trophic_Level=Trophic_Level,IUCN_status=IUCN_status,Temp_min=temp_min,Temp_max=temp_max, Repro.Mode=Repro.Mode, Repro.Fertil=Repro.Fertil, Repro.MatingT=Repro.MatingT, Repro.SpawnFreq=Repro.SpawnFreq, Repro.ParentCare=Repro.ParentCare, Egg.Shape=Egg.Shape, Egg.Attrib=Egg.Attrib, Egg.Color=Egg.Color, Larvae.LenghtAtBirth, Larvae.PlaceDevelop, Larvae.YolkSac))
    rownames(res)=x  
  }
  
  return(res)
  
}

###  
get_item <- function (x="Regalecus-glesne",page="Synonyms",server=".se"){
  
  url2 <- paste("https://www.fishbase",server,"/summary/",x,".html",sep="")
  c <- htmlParse(getURLContent(url2,followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE )))
  link_list <- getHTMLLinks(c,externalOnly=TRUE,xpQuery="//a/@href", baseURL = docName(c))
  
  if(length(link_list)==0){
    return("Not in Fishbase")
  } else {
    
    
    
    if(page=="Morphology"){
      if(length( grep (pattern = "MorphDataList",link_list))==0){
        return("No Morphology information for this species")
      } else {
        link_url <-  substr(link_list [[grep (pattern = "MorphDataList",link_list)[1]]],4,nchar(link_list [[grep (pattern = "MorphDataList",link_list)[1]]]))
        #link_url <- gsub("MorphDataList","MorphDataSummary",link_url)
        #link_url <- gsub("GenusName","genusname",link_url)
        #link_url <- gsub("SpeciesName","speciesname",link_url)
        #a <- strsplit(link_url,"?",fixed=T) 
        #b <- strsplit(a[[1]][2],"&",fixed=T)[[1]][2:3]
        #link_url <- paste0(a[[1]][1],"?",paste(b[1],b[2],sep="&"),"&autoctr=18")
        
      } # end of if
    }
    
    
    
    
    if(page=="Age/Size"){  
      if(length( grep (pattern = "PopCharList",link_list))==0){
        return("No Age/Size information for this species")
      } else {
        link_url <-  substr(link_list [[grep (pattern = "PopCharList",link_list)[1]]],2,nchar(link_list [[grep (pattern = "PopCharList",link_list)[1]]]))
      } # end of if
    }
    
    if( page=="Synonyms"){
      link_url <-  substr(link_list [[grep (pattern = "Synonyms",link_list)[1]]],4,nchar(link_list [[grep (pattern = "Synonyms",link_list)[1]]]))
    } # end of if
    
    if(page=="Reproduction"){
      
      if(length( grep (pattern = "FishReproSummary",link_list))==0){
        return("No reproduction information for this species")
      } else {
        link_url <-  substr(link_list [[grep (pattern = "FishReproSummary",link_list)[1]]],4,nchar(link_list [[grep (pattern = "FishReproSummary",link_list)[1]]])) 
      } # end of ifelse
      
    }
    
    if(page=="Morpho"){
      
      if(length (grep (pattern ="MorphMetList", link_list))==0){
        return("No morpho information for this species")
        
      } else{
        link_int <-  substr(link_list [[grep(pattern = "MorphMetList",link_list)[1]]],4,nchar(link_list [[grep(pattern = "MorphMetList",link_list)[1]]])) 
        urlint <- paste ("https://www.fishbase.se/",link_int,sep="")
        
        pageint <- htmlParse(getURLContent(urlint, followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE )))
        link_listint <- getHTMLLinks(pageint, externalOnly = TRUE, xpQuery = "//a/@href", baseURL = docName(c))
        
        
        if(length (grep (pattern ="picname", link_listint))==0){
          return("No pictures for this species")
        } else {
          link_url <-  link_listint[grep("picname",link_listint)]    
        }
        
      }
      
    } # end of if
    
    if (page=="LW"){
      link_url <-  substr(link_list [[grep (pattern = "LW",link_list)]],2,nchar(link_list [[grep (pattern = "LW",link_list)]]))  	
    } # end of if
    
    if( page=="FoodItems"){
      link_url <-  substr(link_list [[grep (pattern = "FoodItemsList",link_list)]],4,nchar(link_list [[grep (pattern = "FoodItemsList",link_list)]]))
    } # end of if
    
    url_end <- paste ("https://www.fishbase.se/",link_url,sep="")
    
    if(page=="FoodItems"){
      tt <- getURLContent(url_end, followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
      res <- readHTMLTable(tt,header=NA,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer())
      
      return(res$dataTable)
    } # end of if
    
    if(page=="Reproduction"){
      tt <- getURLContent(url_end, followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
      res <- readHTMLTable(tt,header=FALSE,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer())
      
      return(res$dataTable)
    } # end of if
    
    if(page=="Synonyms"){
      tt <- getURLContent(url_end, followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
      res <- readHTMLTable(tt,header=NA,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer())
      
      res <- res[[1]]
      colnames(res)<- c("Warning","Synonym","Author","CoLStatus","Valid", "Synonymy", "Combination")
      res <-res[-1,]
      return(res)
    } # end of if 
    
    if(page=="LW"){
      tt <- getURLContent(url_end, followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
      res <- readHTMLTable(tt,header=NA,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer())
      
      res <-res[[3]]
      selected_length <- names(which.max(summary(res[,"Lengthtype"])))
      res <- res[which(res$Lengthtype==TL),]
      sup_lines <- which(res[,4]=="yes")
      
      if(length(sup_lines)!=0){
        res <- res[-sup_lines,]
      } # end of if
      
      a <- mean(as.numeric(as.character(res[,"a"])))
      b <- mean(as.numeric(as.character(res[,"b"])))
      Length_type <- selected_length
      
      n<-as.numeric(as.character(res$n))
      n<-sum(n[!is.na(n)])
      
      return(data.frame(a=a,b=b,Length_type=Length_type,n=n))
    } # end of if
    
    
    if(page=="Morphology"){
      tt <-htmlParse(readLines(base::url(url_end)))
      res <- readHTMLTable(tt,header=FALSE,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer())
      return(res[[3]][-9,])
      
    }
    
    
    if(page=="Age/Size"){
      tt <- getURLContent(url_end, followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
      res <- readHTMLTable(tt,header=NA,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer())
      return(res)
    }
    
    
    
    #res<-yopyep[[1]]
    if(page=="Morpho"){
      if(length(url_end)>1){
        
        yop <- lapply(1:length(url_end), function(v){getURLContent(url_end[v], followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))})
        res <- do.call(cbind,lapply(yop,function(z){readHTMLTable(z,header=NA,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer())[[1]]}))
        
      } else {
        tt <- getURLContent(url_end, followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE ))
        res <- readHTMLTable(tt,header=NA,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer())[[1]]
        
        #Aspect_ratio_of_caudal_fin is not in pixels nor cm
      } 
      
      return(res)
      
    }
    
  } # end of if else
  
}


### to format Morpho tables
### 
morphomet <- function(x=all_fam_morpho[[65]],sp_names=names(all_fam_morpho[65])) {
  if(dim(x)[2]==2) {
    
    res <- x[-c(1,2,3,16),]
    rownames(res) <- gsub(" ","_",res[,1])
    
    Fsize <- do.call(rbind,strsplit(as.character(res[1,2]),split=" "))[,1]
    F<- do.call(rbind,strsplit(as.character(res[2:11,2]),split=" "))[,1]
    Traits <- c(Fsize,F,as.numeric(as.character(res[12,2])))
    
    Traits<-as.numeric(as.character(Traits))
    
    Traits[2:9] <- (Traits[1]/100)*Traits[2:9]
    Traits[10:11] <- (Traits[9]/100)*Traits[10:11]
    Traits[13] <- as.character(x[2,2])
    
    ret <- data.frame(t(Traits),sp_names,colnames(x)[2],stringsAsFactors=FALSE)
    
    colnames(ret)<- c(rownames(res),"Sex","genus_species","picture_ID")
    rownames(ret)<-paste(sp_names,colnames(res)[2],sep="")
    
    ret
    
  }else{
    
    rownames(x) <- gsub(" ","_",x[,1])
    Columns<-grep("Picture Used",colnames(x))
    
    res_3<-x[,-c(grep("Picture Used",colnames(x)))]
    
    res_4 <- res_3[-c(1,2,3,16),]
    
    Fsize_pixel <- apply(res_4[1,],1, function(y) {do.call(rbind,strsplit(as.character(y),split=" "))})[1:(dim(x)[2]/2),1]
    
    if(length(grep(",",Fsize_pixel))>0){Fsize_pixel <- gsub(",","",Fsize_pixel)}
    
    Fsize_SL <- apply(res_4[2,],1, function(y){do.call(rbind,strsplit(as.character(y),split=" "))})[1:(dim(x)[2]/2),1]
    
    Ftraits <-apply(res_4[3:12,],2, function(z){do.call(rbind,strsplit(as.character(z),split=" "))[,1]} )
    if(nrow(Ftraits)!=nrow(res_4[3:12,])){Ftraits <- rbind(Ftraits,c(NA,NA))}
    
    Traits <- rbind (Fsize_pixel,Fsize_SL,Ftraits)
    
    Traits<-apply(Traits,2, function (j) {as.numeric(as.character(j))})
    
    Traitsfinal<-apply(Traits,2,pixel_calculation)
    
    MeanFtraits <- apply(Traitsfinal,1,function(w){ mean(w,na.rm=TRUE) })
    
    res <- data.frame(rbind(t(Traitsfinal),MeanFtraits))
    colnames(res)<-rownames(res_4)
    
    rownames(res)[dim(res)[1]] <- paste(sp_names,"MeanFtraits",strsplit(rownames(res)[2],"_")[[1]][1],sep="_")
    
    rownames(res)[1:(dim(res)[1]-1)]<-paste(sp_names,rownames(res)[1:(dim(res)[1]-1)],sep="_")
    res<- cbind(res,Sex=c(t(res_3[2,]),"NA"),genus_species=c(rep(sp_names,dim(res)[1])),picture_ID=rownames(res))
    
    res
    
  } #end of else
  
} #end of function 

pixel_calculation <- function (x=Traits[,1]) {
  
  x[2:9] <- (x[1]/100)*x[2:9]
  x[10:11] <- (x[9]/100)*x[10:11] 
  
  x
} #end of function to calculate pixels from % of TL and HL


###extract names of all species per family from FB:
###
get_sp_familyFB <- function(x="Labridae"){
  
  url2<-paste ("https://www.fishbase.org/Nomenclature/FamilySearchList.php?Family=", x,sep="") #construction du chemin url
  tt <- getURLContent(url2, followlocation=TRUE,.opts =  RCurl::curlOptions(ssl.verifypeer=FALSE )) #lecture de la page
  data.frame(readHTMLTable(tt,header=NA,colClasses=NULL,skip.rows=integer(),trim =TRUE,elFun=xmlValue,as.data.frame=TRUE,which=integer()))
  
}

sp_per_family <- as.character(get_sp_familyFB(x="Sparidae")[,1])
sp_per_family<-gsub(" ","-",sp_per_family)
