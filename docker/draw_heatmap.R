

# Col.Groups need to be a list of groups
Draw.Heatmap<-function(Sig.nData=NA,  Col.Groups=NA,  Row.Groups=NA,
                       title="", scale="none", dendrogram="both",
                       col.grad=50,
                       col.panel=c("blue", "white", "red"),
                       outfile=NA,
                       margin.x=20, margin.y=10, Rowv=TRUE,
                       width=2400, height=2400, res=150){
  r.colors=c("grey", "green", "yellow", "blue", "red")
  c.colors=c("grey", "blue", "red", "green", "cyan")
  Col.Group=if(is.na(Col.Groups)){colnames(Sig.nData)
  }else{as.character(unlist(Col.Groups))}
  Row.Group=if(is.na(Row.Groups)){rownames(Sig.nData)
  }else{as.character(unlist(Row.Groups))}
  #############
  #    Sig.nData=Probe.w.Signal.mtx
  CurrMtx=Sig.nData[Row.Group, Col.Group ]



  # Assign Sample colors

  sample.ls=colnames(CurrMtx)
  samplecolors=sample.ls
  if(!is.na(Col.Groups)){
    for(i in 1:length(names(Col.Groups))){
      Curr.Group=names(Col.Groups)[i]
      curr.entr.group=Col.Groups[[Curr.Group]]
      samplecolors[samplecolors %in% curr.entr.group]=c.colors[(i %% length(c.colors))+1]
    }
  }else{
    samplecolors[1:length(samplecolors)]="blue"
  }

  row.ls=rownames(CurrMtx)
  rowcolors=row.ls
  if(!is.na(Row.Groups)){
    for(i in 1:length(names(Row.Groups))){
      Curr.Group=names(Row.Groups)[i]
      curr.entr.group=Row.Groups[[Curr.Group]]
      rowcolors[rowcolors %in% curr.entr.group]=r.colors[(i %% length(r.colors))+1]
    }
  }else{
    rowcolors[1:length(rowcolors)]="green"
  }


  if(!(is.na(outfile))){
    png(file=outfile, width=width, height=height, res=res)
  }
  # colorpanel(10, "green", "white", "red")
  heatmap.2(CurrMtx, col=colorpanel(col.grad, col.panel[1],
                                    col.panel[2], col.panel[3]), trace="none",
            main=title,
            #ColSideColors=samplecolors,
            #RowSideColors=probecolors,
            xlab="Sample",
            ylab="",
            dendrogram=dendrogram,
            scale=scale,
            #scale="column",
            #scale="row",
            margin=c(margin.x, margin.y),
            Rowv=Rowv,
            Colv=if(dendrogram=="row"){FALSE}else{TRUE})
  if(!(is.na(outfile))){
    dev.off()
  }
}
