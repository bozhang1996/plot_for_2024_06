setwd("D:/submission/others/energy/z25/Figure6/netcir")
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(reshape2)
library(ggplot2)
library(psych)
data<- read.csv("data.csv",header = TRUE,row.names = 1,check.names=FALSE,comment.char = "")

type = read.csv("type.csv",header = TRUE,check.names = FALSE)
head(type)

cor <- psych::corr.test(data, use = "pairwise",
                        method="spearman",
                        adjust="fdr", 
                        alpha=0.10,
                        ci=T) # ci=FALSE,不进行置信区间计算，数据量较大时，可以加快计算速度。
cor.r <- data.frame(cor$r) # 提取R值
cor.p <- data.frame(cor$p) # 提取p值
colnames(cor.r) = rownames(cor.r)
colnames(cor.p) = rownames(cor.p) # 变量名称中存在特殊字符，为了防止矩阵行名与列名不一致，必须运行此代码。

# write.csv(cor.r,"cor.r.csv",quote = FALSE,col.names = NA,row.names = TRUE) # 保存结果到本地
# write.csv(cor.p,"cor.p.csv",quote = FALSE,col.names = NA,row.names = TRUE)

head(cor.r)
head(cor.p)

#-----保留p值及计算相关性表格#----
cor.r$from = rownames(cor.r) 
cor.p$from = rownames(cor.p)
p = cor.p %>% 
  gather(key = "to", value = "p", -from) %>%
  data.frame()

cor.data = cor.r %>% 
  gather(key = "to", value = "r", -from) %>%
  data.frame() %>%
  left_join(p, by=c("from","to")) %>%
  filter(p <= 0.05, from != to) %>%
  mutate(
    linecolor = ifelse(r > 0,"positive","negative"), # 设置链接线属性，可用于设置线型和颜色。
    linesize = abs(r) # 设置链接线宽度。
  ) # 此输出仍有重复链接，后面需进一步去除。
head(cor.data)



# 3.1 准备节点属性文件
## 3.1.1 计算每个节点具有的链接数(degree)
c(as.character(cor.data$from),as.character(cor.data$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("name", "n")

## 3.1.2 添加变量分类属性
vertices <- vertices %>%
  select(-n) %>% # 因为此处的n不准确，所以先删除。
  left_join(type,by="name")

## 3.1.3 对节点属性表进行排序
#网络图中节点会按照节点属性文件的顺序依次绘制，
##为了使同类型变量位置靠近，按照节点属性对节点进行排序。
##环境变量最好置于表格末尾，否则后续绘图可能会效果不好
##此处使用type进行排序，依次是细菌、真菌和环境因子,
##然后按照group排序，使同一属的物种彼此靠近。
vertices$type <- factor(vertices$type,
                        levels = c("Firmicutes","Bacteroidetes","Actinobacteria","Deferribacteres","Verrucomicrobia"
                        ,"Proteobacteria","TRP","5HT","Indole" ,"Kynurenine" ))
vertices <- vertices %>%
  arrange(type,group) 

dim(vertices)
head(vertices)


## 3.2.1 构建graph数据结构
graph <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
graph # 非简单图

## 3.2.2 转为简单图，去除cor.data数据表中的重复链接
is.simple(graph) # 非简单图，链接数会偏高，所以需要转换为简单图。
E(graph)$weight <- 1 # 将链接权重赋值为1
graph <- igraph::simplify(graph,edge.attr.comb = "first")# 转为简单图
is.simple(graph)
E(graph)$weight <- 1 # 将链接权重赋值为1
is.weighted(graph)
graph
V(graph)$degree <- degree(graph)
net.data  <- igraph::as_data_frame(graph, what = "both")$edges # 提取链接属性
vertices  <- igraph::as_data_frame(graph, what = "both")$vertices # 提取节点属性


unique(vertices$group) # 7个分组信息
#??ggsci # 查看ggsci包帮助文件，选择颜色。
library(ggsci)
library(scales)
mycolor = pal_d3("category20",alpha = 1)(20)
mycolor
cols = mycolor[c(1:7,9:13,14)]
cols <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
          "#4aef7b", "#e86502", "#9ed84e", "#8249aa", "#99db27", "#e07233", "#ff523f",
          "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
          "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
          "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
          "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
show_col(cols) # 以图形形式展示颜色

## 链接线颜色
#linetypes = c("positive" ="solid","negative" ="dashed") # 设置线类型
color = c("positive" ="#D62728FF","negative" ="#2CA02CFF")

#-----生成角度并构建graph#-----
bar1 = nrow(vertices[vertices$group == "metabolites",])
bar1 # 环境因子
bar2 = nrow(vertices[vertices$group != "metabolites",])
bar2 # 物种


vertices$id =1
vertices$id[vertices$group == "metabolites"] = seq(1, bar1)
vertices$id[vertices$group != "metabolites"] = seq(1, bar2)
head(vertices)


vertices$var_angel <- NA
vertices$var_angel[vertices$group == "metabolites"] <- 360 * (vertices$id[vertices$group == "metabolites"]-0.5)/bar1 
vertices$var_angel[vertices$group != "metabolites"] <- 360 * (vertices$id[vertices$group != "metabolites"]-0.5)/bar2


vertices$x <- ifelse(vertices$group == "metabolites",cos(vertices$var_angel)/2,cos(vertices$var_angel))
vertices$y <- ifelse(vertices$group == "metabolites",sin(vertices$var_angel)/2,sin(vertices$var_angel))
head(vertices)

vertices$hjust <- ifelse(vertices$var_angel > 180, 1, 0)
vertices$angle <- ifelse(vertices$var_angel > 180, 90-vertices$var_angel+180, 90-vertices$var_angel)

vertices$hjust[vertices$group != "metabolites"] = "center" # 或0.5
vertices$angle[vertices$group != "metabolites"] = 0

net.data <- net.data %>% 
  left_join(vertices,by=c('from'='name')) 

graph <- graph_from_data_frame(net.data, vertices = vertices, directed = FALSE)
graph

#write.csv(vertices,"relation.csv", sep = "\t")
write.csv(net.data,"netdata0911.csv",quote = FALSE,col.names = NA,row.names = FALSE) # 保存结果到本地。
layout <- create_layout(graph,layout = 'linear', circular = TRUE)
head(layout)
layout$x <- vertices$x
layout$y <- vertices$y
head(layout)

net.cir1 <- ggraph(layout) +
  geom_edge_link(aes(edge_colour = as.factor(linecolor),
                     edge_width = abs(r)
                     #edge_linetype = linetype # 想用不同线形表示相关性正负，取消此句注释
  ),
  #edge_width=0.5,
  edge_alpha=0.9 ) + # 使用直线作为链接线。
  ###theme()中的legend参数不利于多个图例设置，这里在scale*函数中直接设置图例参数###
  scale_edge_colour_manual(values=color,
                           breaks = c("positive","negative"),
                           guide = guide_legend(title="Correlation",
                                                direction="vertical",
                                                order=3,# 图例排列序号
                                                ncol=1,
                                                byrow=FALSE,
                                                title.theme = element_text(
                                                  size = 14,
                                                  face = "bold",
                                                  colour = "black"),
                                                label.theme = element_text(
                                                  size = 12,colour = "black")
                           )) +
  scale_edge_width(
    breaks = seq(0.2,1,0.2),
    label = seq(0.2,1,0.2),
    range = c(0.2,1),
    guide = guide_legend(title="Pearson |r|",
                         direction="vertical",
                         order=1,# 图例排列序号
                         ncol=1,
                         byrow=FALSE,
                         title.theme = element_text(
                           size = 14,
                           face = "bold",
                           colour = "black"),
                         label.theme = element_text(
                           size = 12,colour = "black")
    ))+
  ###想用不同线形表示相关性正负，取消此句注释###
  #scale_edge_linetype_manual(values = linetypes,
  #name = "Correlation",
  #breaks = c("positive","negative"))+
  geom_node_point(aes(size=degree, fill=as.factor(group)),
                  shape=21,alpha=1) +
  #scale_size_continuous(name = "Degree",range=c(10,15)) +
  scale_size(
    breaks = seq(3,max(vertices$degree),3),
    label = seq(3,max(vertices$degree),3),
    range = c(5, 15),
    guide = guide_legend(title="Degree",
                         direction="vertical",
                         order=3,# 图例排列序号
                         ncol=2, # 图例排两列
                         byrow=FALSE,
                         title.theme = element_text(
                           size = 14,
                           face = "bold",
                           colour = "black"),
                         label.theme = element_text(
                           size = 12,colour = "black")
    ))+
  ###可使用cols = c(...,"env"="black")为每个节点设置指定颜色，适合节点分类少的网络图###
  scale_fill_manual(values=cols,
                    guide = guide_legend(title="Genus",
                                         #keywidth=grid::unit(2,"cm"),
                                         #keyheight=grid::unit(2,"cm"),
                                         direction="vertical",
                                         order=4,# 图例排列序号
                                         ncol=2,
                                         byrow=FALSE,
                                         title.theme = element_text(
                                           size = 14,
                                           face = "bold",
                                           colour = "black"),
                                         label.theme = element_text(
                                           size = 12,colour = "black",face = "italic")
                    )) +
  ###展示物种+环境因子节点标签，不展示节点标签，则注释掉geom_node_text()函数###
  #geom_node_text(aes(x = x*1.1, y=y*1.1,
  #label=as.character(name), 
  #angle=angle,hjust=hjust,
  #color=as.factor(group)),size=4.5,
  #show.legend = FALSE) +
  ###只展示环境因子节点标签###
  geom_node_text(aes(x = x, y=y,
                     label=ifelse(V(graph)$group == "metabolites",as.character(name),""), 
                     angle=angle,hjust=hjust,
                     color = as.factor(group)
  ),
  size=4.75,# 14pt
  show.legend = FALSE) +
  ###将env分类对应的颜色改为“black”,更改节点标签颜色###
  scale_color_manual(values = c(cols[1:6],"black"))+ 
  #expand_limits(x = c(-1, 1), y = c(-1, 1))+ # 设置网络图区域坐标范围
  coord_fixed()+
  theme_minimal()+
  theme(
    legend.title.align = 0,# 图例标题左对齐方式(0,左对齐；1，右对齐；0.5，居中对齐)
    legend.text.align = 0,# 图例文字左侧对齐方式(0,左对齐；1，右对齐；0.5，居中对齐)
    ###图例的左上角位于图中(1,0.9)坐标处###
    #legend.position = c(1,0.9), # 坐标设置图例位置左→右：0-1；下→上：0-1.
    #legend.justification = c(0,1), # 更改图例锚点(左→右：0-1；下→上：0-1,默认(0,1))
    legend.spacing.x = unit(0.1,"cm"), # 图例中各个元素(key与key标签)的距离
    legend.spacing.y = unit(0.1,"cm"),
    legend.box.margin = margin(t=0.2,r=0.2,b=0.2,l=0.2,"cm"), # 设置图例中心与图例框边缘的距离
    plot.margin=unit(rep(0,4), "null"), # 设置图片区与边界的距离,值越大图与图片边界距离越远。
    panel.spacing=unit(c(0,0,0,0), "null"),
    ###去除图中坐标轴和网格线###
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks =element_blank(),
    axis.text =element_blank(),
    axis.title = element_blank()
  )
net.cir1


layout1 <- create_layout(
  induced.subgraph(graph,c(V(graph)$name %in% vertices$name[vertices$group != "metabolites"])), 
  layout = 'linear', circular = TRUE)
layout1

### 计算环境因子子网络节点布局数据表
layout2 <- create_layout(
  induced.subgraph(graph,c(V(graph)$name %in% vertices$name[vertices$group == "metabolites"])), 
  layout = 'linear', circular = TRUE)
layout2[1:2] <- layout2[1:2]/2 # 内环变量的坐标同时除2。


## 4.3.2 create_layout()创建布局数据，并将其中的(x,y)替换为上面计算的坐标。
layout <- create_layout(graph,layout = 'linear', circular = T)
head(layout)

## 4.3.3 替换坐标-要求name的顺序一致。
###rbind(layout2,layout1)，不能用于设置ggraph(graph,layout = layout),绘图时会缺少物种链接线布局信息。
layout$x <- rbind(layout1,layout2)$x
layout$y <- rbind(layout1,layout2)$y
head(layout)

Correlation <- guide_legend(title="Correlation",
                            direction="vertical",
                            order=2,# 图例排列序号
                            ncol=1,
                            byrow=FALSE,
                            title.theme = element_text(
                              size = 14,
                              face = "bold",
                              colour = "black"),
                            label.theme = element_text(
                              size = 12,colour = "black")
)

width_legend <- guide_legend(title="Pearson |r|",
                             direction="vertical",
                             order=1,# 图例排列序号
                             ncol=1,
                             byrow=FALSE,
                             title.theme = element_text(
                               size = 14,
                               face = "bold",
                               colour = "black"),
                             label.theme = element_text(
                               size = 12,colour = "black")
)

size_legend <- guide_legend(title="Degree",
                            direction="vertical",
                            order=3,# 图例排列序号
                            ncol=2, # 图例排两列
                            byrow=FALSE,
                            title.theme = element_text(
                              size = 14,
                              face = "bold",
                              colour = "black"),
                            label.theme = element_text(
                              size = 12,colour = "black")
)

fill_legend <- guide_legend(title="Genus",
                            direction="vertical",
                            order=4,# 图例排列序号
                            ncol=2,
                            byrow=FALSE,
                            title.theme = element_text(
                              size = 14,
                              face = "bold",
                              colour = "black"),
                            label.theme = element_text(
                              size = 12,
                              colour = "black",face = "italic")
)

set_graph_style(plot_margin = margin(0,0,1,0))

net.cir2 <- ggraph(layout) +
  geom_edge_link(aes(edge_colour = as.factor(linecolor),
                    edge_width = abs(r)
  ),
  edge_alpha=0.6 ) + # 改用曲线链接。
  scale_edge_colour_manual(values=color,
                           breaks = c("positive","negative"),
                           guide = Correlation) +
  scale_edge_width(
    breaks = seq(0.2,1,0.2),
    label = seq(0.2,1,0.2),
    range = c(0.2,1),
    guide = width_legend )+
  geom_node_point(aes(size=degree, fill=as.factor(type)),
                  shape=21,alpha=1) +
  scale_fill_manual(values=cols,
                    guide = fill_legend) +
  scale_size(
    breaks = seq(3,max(vertices$degree),3),
    label = seq(3,max(vertices$degree),3),
    range = c(5, 15),
    guide = size_legend)+
  geom_node_text(aes(x = x, y=y,
                     label= ifelse(V(graph)$group == "metabolites",as.character(name),""),
                     angle=angle,hjust=hjust,
                     #color = as.factor(group)
  ),
  color = "black",# 所有节点标签设置为黑色。
  size=4.75,# 14pt
  show.legend = FALSE) +
  #添加物种标签，取消下面代码的注释即可##
  geom_node_text(aes(x = x*1.4, y=y*1.4, # 物种与环境因子标签的放置方式不同，最好分别设置，
  label= ifelse(V(graph)$group != "metabolites",as.character(name),""),
  angle=angle,hjust=hjust,
      ),
  color = "black",# 所有节点标签设置为黑色。
  size=4.75,# 14pt
  show.legend = FALSE) +
  scale_color_manual(values = c( "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92")) +
  coord_fixed()+
  theme_graph()
net.cir2

cairo_pdf("netcir2.pdf",height=12,width=12,family="Times")
print(net.cir2) 
dev.off()
