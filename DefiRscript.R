####Importation des données####

#Pour utiliser ce script, une table de comptage doit etre utilise

expData <- read.table("Mito_Genes.txt", row.names = 1, sep = "\t", header = T)
if(exists("expData")==FALSE){
  print("Veuillez charger une table de comptage")
}

expData2 <- read.table("Mito_Genes2.txt", row.names = 1, sep = "\t", header = T)
#Matrice avec des NA

####Fonction Gaëlle représentation expression gene####
plotGenes <- function(expData, title = "", yMin = 0, yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(floor(yMin), ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}


#####Condition de la méthode de distance####

#Pour realiser un clustering, il est preferable de creer des matrices de distances. En utilisant une table de comptage, deux choix sont proposés, un calcul de distance euclidienne ou de correlation. 

#Si l'operateur a bien choisi les parametres une matrice de distance est enregistre dans son environnement. Sinon un message specifie les choix possibles

#Creation fonction
Methode_Distance<-function(expData,MetDist="euclidienne"){
  #Controle des parametres d'entree  
  if(MetDist=="euclidienne"){
    #utilisation de la fonction dist pour creer une matrice de distance euclidienne  
    matDist <- dist(expData, method = "euclidean")
    
  }else if(MetDist=="correlation"){
    #cor() permet de creer une matrice de correlation pour la convertir en matrice de distance, les valeurs doivent etre positives et le tableau transpose t()  
    matDist <- as.dist(1 - cor(t(expData)))
    
  }else{
    print("Veuillez choisir une méthode de distance : euclidienne ou correlation")
  }
  assign("matDist",matDist,envir = parent.frame()) 
}


####Condition de la méthode de clustering####
#Une fois la matrice de distance realise, le clustering pourrait etre realise. Le parametre d'entree doit etre une matrice de distance (meme si la fonction kmeans peut utiliser directement une table de comptage mais qui sera transforme en matrice de distance euclidienne). Pour cette fonction, le nombre de cluster souhaite doit egalement etre renseigne. Deux choix de methode sont possibles, kmeans (parametre par defaut) et clustering hierarchique (HCL). Pour HCL, plusieurs methodes sont possibles (voir hclust()) avec ward.D2 par defaut.

#Si l'operateur a bien choisi les parametres le clustering est enregistre dans son environnement. Sinon un message specifie les choix possibles 

#Creation fonction
Methode_Clustering<-function(matDist,N,MetClust="kmeans",MetClustHier="ward.D2"){
  #Controle des parametres d'entree 
  if(MetClust=="kmeans"){
    #utilisation de la fonction kmeans  pour un clustering de "k-moyens"   
    resClust <- kmeans(matDist, centers = N)
    
  }else if(MetClust=="HCL"){
    #utilisation de la fonction hclust() pour le clustering hierarchique  
    resClust <- hclust(matDist, method = MetClustHier)
    
  }else{
    print("Veuillez choisir une méthode de clustering : kmeans ou HCL")
  }
  assign("resClust",resClust,envir = parent.frame()) 
}

####Création du vecteur de cluster pour le graphique####

#En fonction du choix de la methode de clustering, les sorties de fonctions sont differentes. L'assignation des genes et de leur cluster doit etre stocke dans un vecteur pour pouvoir realiser un graphique. 
#Pour la methode kmeans, la sortie est un objet R et les vecteurs sont contenus dans $cluster.
#Pour la methode HCL, la sortie est un dendogramme. La fonction cutree permet de couper à partir d'un nombre de ramification donnee et donc de cluster.

#Creation fonction
vector_cluster<-function(resClust,N,MetClust="kmeans"){
  #Controle des parametres d'entree 
  if(MetClust=="kmeans"){
    #recuperation de l'objet cluster  
    Clust <- resClust$cluster
    
  }else if(MetClust=="HCL"){
    #utilisation de cutree() pour couper le dendogramme et recupere le nombre de cluster souhaite et leurs genes assignes  
    Clust <- cutree(resClust,k=N) 
    
  }
  assign("Clust",Clust,envir = parent.frame()) 
}  


####Condition du graphique####

#La fonction graphique va generer deux type de graphique "profil d'expression" (profil) ou heatmap avec la possibilite d'afficher les deux (Both).
#Si ils sont affiche separemment, ces graphiques seront generes par des fonctions classique de R.
#Si les deux graphiques sont affiches, ggplot sera utilise avec la possibilite de sauvegarde du graphique en PNG

#Creation fonction
Graphique<-function(Clust,N,MetGraph="profil",ggsave=FALSE,width=15,height=10){
  library(RColorBrewer)
  #Librarie contenant des couleurs plus sympa pour la heatmap
  library(dplyr)
  #librarie pour manipuler des donnees tidy
  library(tibble)
  #librairie pour modifier des tableau tidy
  library(tidyr)
  #librairie pour mettre en "ordre" sur ggplot
  library(ggplot2)
  #librairie de graphique
  library(gridExtra)
  #librairie pour modifier les grilles de graphiques et d'en contenir sur la meme fenetre
  
  #coul<-colorRampPalette(brewer.pal(8,"PiYG"))(25)
  
  #Controle des parametres d'entree 
  if(MetGraph=="profil"){
    
    #boucle pour generer un graphique par cluster  
    for(i in 1:N){
      cluster <- expData[which(Clust == i),]
      #Utilisation de la fonction graphique de Gaelle
      plotGenes(cluster, yMax = 100,title=paste("Profil d'expression Cluster",i))
    }
    
    
  }else if(MetGraph=="heatmap"){
    
    for(i in 1:N){
      cluster <- expData[which(Clust == i),]
      #utilisation de la fonction heatmap() avec une palette de couleur
      heatmap(as.matrix(cluster),main=paste("Heatmap Cluster",i),col=colorRampPalette(brewer.pal(8,"RdBu"))(25))
    }
    
  }else if(MetGraph=="Both"){
    
    for(i in 1:N){
      cluster <- expData[which(Clust == i),]
      tidydata<-cluster %>%
        #Donner un nom à la colonne contenant les noms de lignes
        rownames_to_column() %>%
        #transformer le tableau en tidy a partir de la colonne 2
        pivot_longer(c(2:ncol(.)))%>%
        #grouper les noms de gene entre eux
        group_by(name) %>%
        #creer une nouvelle colonne avec un numero pour chaque nom de gene (necessaire pour certains types de                graphique avec des variables continus)
        mutate(group_id = cur_group_id())
      
      #Creation d'une variable avec le graphique de type profil d'expressions car en cas d'utilisation de ggsave(), l'aperçu du graphique n'apparait pas
      #geom_line permet de voir les profils d'expressions
      Profil<-ggplot(tidydata,aes(group_id,value))+geom_line(aes(group=rowname),color="grey")+
        #stat_summary affiche la moyenne des profls
        stat_summary(aes(y = value,group=1), fun=mean, colour="red", geom="line",group=1,linetype="dashed")+
        #le theme doit être modifier en retirant les valeurs de certains axes pour plus de lisibilté
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        #Ajout du titre et nom des axes
        xlab("Time point")+ ylab("Gene expression level")+  ggtitle(paste("Expression Profile Cluster ",i))+
        #Permet aux axes du graphiques de rester dans des valeurs positives
        scale_x_continuous(expand = c(0, 0), limits = c(NA,NA))
      
      #Creation d'une variable avec le graphique de type heatmap car en cas d'utilisation de ggsave(), l'aperçu du graphique n'apparait pas
      #geom_tile permet une heatmap
      Heat<-ggplot(tidydata,aes(name,rowname,fill=value))+geom_tile()+
        #Ajout du nom des axes  
        xlab("Time point")+ ylab("Gene expression level")+
        #le theme doit être modifier en retirant les valeurs de certains axes pour plus de lisibilté
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
        #Evite l'overlapping du noms des variables sur les axes.  
        scale_y_discrete(guide = guide_axis(check.overlap = TRUE))
      
      
      Graph<-grid.arrange(Profil,Heat,ncol=2,nrow=1)
      
      if(ggsave==TRUE){
        
        ggsave(paste("Cluster",i,".png"),Graph,width=width,height=height)
        
      }
    }
  }else{
    print("Veuillez choisir un type de graphique : profil ou heatmap")
  }
}


##### Fonction Main the big ONE #####

#la fonction BigOne va permettre de faire appel a toutes les fonctions precedemments en partant d'une matrice de comptage et un nombre de cluster 

BigOne<-function(expData,N,MetDist="euclidienne",MetClust="kmeans",MetGraph="profil",MetClustHier="ward.D2",ggsave=FALSE){
  Methode_Distance(expData,MetDist=MetDist)
  Methode_Clustering(matDist,N=N,MetClust=MetClust,MetClustHier=MetClustHier)
  vector_cluster(resClust,N=N,MetClust=MetClust)
  Graphique(Clust,N=N,MetGraph=MetGraph,ggsave=ggsave)
}

##### Fonction Main the big ONE v2 #####

#la fonction v2 va faire appel a BigOne mais en ajouter une option de correction si des NA sont presents dans la matrice.
#En cas de NA, un menu sera affiche pour demander si l'utilisateur veut corriger sa matrice.
#Si une ligne comprend plus de 80% de NA, elle sera supprimee et l'utilisateur sera informé du numero de ligne correspondant. 
#Si la ligne comprend moins de 80%, la methode des plus proches voisins sera utilisee par la fonction knn() du package VIM.
#la fonction knn() semble gourmande et donc desactive dans ce script (fonctionne sur les PC du DU).

BigOnev2<-function(expData,N,MetDist="euclidienne",MetClust="kmeans",MetGraph="profil",MetClustHier="ward.D2",ggsave=FALSE){  
  SumNA<-sum(is.na(expData))
  
  if (SumNA>1){
    #la fonction menu() va demander a l'utilisateur si il souhaite corriger sa matrice. L'option graphics va permettre d'aificher une boite de dialoue dans une nouvelle fenetre Windows.
    Correction<-menu(c("Oui", "Non","La réponse D"), graphics=FALSE,title=paste("NA", "détecté(s) dans la matrice, voulez-vous la corriger ?"))
    if(Correction==1){
      expData80percentNA<-nrow(expData[rowSums((is.na(expData)))/ncol(expData) >= 0.8,])
      if(expData80percentNA>=1){
        print("les lignes suivantes contennant plus de 80% de NA vont être supprimées")
        print(paste(which((rowSums((is.na(expData)))/ncol(expData)) >= 0.8) ))
        expData<-expData[rowSums((is.na(expData)))/ncol(expData) <= 0.8,]
        
      }
      
      library(VIM)
      # imputation grâce la méthode des kNN :
      expData<-kNN(expData)
      print("TABERNACLE! La matrice est corrigée, la fonction BigOne va débuter.")
      BigOne(expData,N,MetDist=MetDist,MetClust=MetClust,MetGraph=MetGraph,MetClustHier=MetClustHier,ggsave=FALSE)
    }
    else if(Correction==3){
      print("Satané Pancake") 
    }
    else{
      print("La matrice doit être corrigé avant de continuer")
    }
  }else{
    BigOne(expData,N,MetDist=MetDist,MetClust=MetClust,MetGraph=MetGraph,MetClustHier=MetClustHier,ggsave=FALSE)
  }
  
}