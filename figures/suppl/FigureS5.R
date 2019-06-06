############COMPARISON OF CLASSIFERS###############################
Classifiers <- c(rep("SVM" , 4) , rep("KNN" , 4) , rep("RF", 4) , rep("LD" , 4), rep("DT" , 4))
Metrics=rep(c("accuracy" , "recall" , "precision", "f1 score") , 5)
Values = c(0.942295918367357, 0.8598866279069756, 0.8737637430242022, 0.8632885672322217, 
           0.9341632653061334, 0.8530658914728683, 0.8500497979782959, 0.8474803971536752, 
           0.9384795918367452, 0.8655639534883712, 0.8588009695873556, 0.8586385965752061, 
           0.9376428571428675, 0.882260658914728, 0.8511795465122552, 0.8626623965269422, 
           0.9390918367347042, 0.8594593023255824, 0.8628965798240885, 0.8573559684227736)

data=data.frame(Classifiers,Metrics, Values)

ggplot(data, aes(fill=Metrics, y=Values, x=Classifiers))+ geom_bar(position="dodge", stat="identity") + coord_cartesian(ylim = c(0.85, 0.95)) + theme(text = element_text(size=15))