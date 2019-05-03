#!/usr/bin/env python
docstring='''
obo2csv.py go-basic.obo
    parse OBO format gene ontogology file "go.obo", 
    write tab-eliminated parsed GO Terms to:
    go.obo.MF.csv,   go.obo.BP.csv,   go.obo.CC.csv
    go.obo.MF.list,  go.obo.BP.list,  go.obo.CC.list
    go.obo.MF.is_a,  go.obo.BP.is_a,  go.obo.CC.is_a
    go.obo.MF.has_a, go.obo.BP.has_a, go.obo.CC.has_a
    go.obo.obsolete
    go.obo.GOaltids.dat
    go.obo.uninformative

    Where F, P, C stands for 3 Aspect/namespace of GO:
    molecular_function, biological_process, cellular_component.

    *.csv contains
    id, name, Aspect, is_a (direct parent)

    *.list contains
    id ! name

    *.is_a contains
    id\tis_a (culmulative)

    *.has_a contains
    id, has_a (culmulative)

    only non-obsolete term will be written
    
    In addition, *.obsolete lists obsolete Terms

    *GOaltids.dat contains:
    alt_id   id

    *.uninformative lists terms that should not be used for direct annotation
'''
import sys,os

GO_namespace_to_Aspect={
    "molecular_function":'F', "MF":'F', "F":'F',
    "biological_process":'P', "BP":'P', "P":'P',
    "cellular_component":'C', "CC":'C', "C":'C',
}

GO_Aspect_to_namespace={
    'F':"molecular_function", 'MF':"molecular_function", 
    'P':"biological_process", 'BP':"biological_process",
    'C':"cellular_component", 'CC':"cellular_component",
    "molecular_function":"molecular_function",
    "biological_process":"biological_process",
    "cellular_component":"cellular_component",
}

# comment that shows a term is uninformative
uninformative_txt="Note that this term is in the subset of terms that should not be used for direct"

class GO_Term:
    # class to store none-obsolete GO Terms
    def __init__(self,Term_txt=''):
        self.id=''
        self.name=''
        self.namespace=''
        self.definition='' # def
        self.comment=''

        self.xref=set()
        #self.synonym=[]
        self.alt_id=set()
        self.is_a=set()
        #self.is_obsolete=False

        #self.obo=Term_txt

        for line in Term_txt.splitlines():
            if   line.startswith("id: "):
                self.id        =line[len("id: "):]
            elif line.startswith("name: "):
                self.name      =line[len("name: "):]
            elif line.startswith("namespace: "):
                self.namespace =line[len("namespace: "):]
            elif line.startswith("def: "):
                self.definition=line[len("def: "):]
            elif line.startswith("comment: "):
                self.comment   =line[len("comment: "):]
            elif line.startswith("xref: "):
                self.xref.add(line[len("xref: "):])
            elif line.startswith("alt_id: "):
                self.alt_id.add(line[len("alt_id: "):])
            elif line.startswith("is_a: "):
                self.is_a.add(line[len("is_a: "):])

    def __str__(self):
        '''print GO Term in CSV format'''
        return '\t'.join([
            self.id, 
            self.name, 
            GO_namespace_to_Aspect[self.namespace], 
            #self.definition,
            #', '.join(self.xref), 
            #self.comment,
            '; '.join(self.is_a),
            ])

    def __repr__(self):
        '''print GO Term in OBO format'''
        obo_list=["[Term]"]
        if self.id:
            obo_list.append("id: "+self.id)
        if self.name:
            obo_list.append("name: "+self.name)
        if self.namespace:
            obo_list.append("namespace: "+self.namespace)
        if self.definition:
            obo_list.append("def: "+self.definition)
        if self.comment:
            obo_list.append("comment: "+self.comment)
        if self.xref:
            obo_list+=["xref: "+xref for xref in self.xref]
        if self.alt_id:
           obo_list+=["alt_id: "+alt_id for alt_id in self.alt_id]
        if self.is_a:
           obo_list+=["is_a: "+is_a for is_a in self.is_a]
        return '\n'+'\n'.join(obo_list)+'\n'
    
    def short(self):
        '''print GO Term in short format'''
        return "%s ! %s"%(self.id,self.name)

class obo(dict):
    '''class to store obo format GO ontogology data
    obo["F"] # molecular_function GO Term
    obo["P"] # biological_process GO Term
    obo["C"] # cellular_component GO Term

    obo["F"]["Term"]        # a list of non-obsolete GO Term
    obo["F"]["is_obsolete"] # a list of obsolete GO Term
    obo["F"]["alt_id"]      # a dict of alternative id, 
                            # key is alt_id, value is primary id
    obo["F"]["is_a"]        # GO hierachy
    obo["F"]["uninformative"] # terms that should not be used for 
                              # direct annotation
    '''
    def __init__(self,obo_txt=''):
        self["F"]={"Term":dict(), "is_a":dict(), "alt_id":dict(),
                   "is_obsolete":[], "uninformative":[]}
        self["P"]={"Term":dict(), "is_a":dict(), "alt_id":dict(),
                   "is_obsolete":[], "uninformative":[]}
        self["C"]={"Term":dict(), "is_a":dict(), "alt_id":dict(),
                   "is_obsolete":[], "uninformative":[]}
        self.append(obo_txt)
        return

    def append(self,obo_txt='',update_hierachy=True):
        '''add new Term to obo class using obo format text "obo_txt"
        update_hierachy: whether to update is_a hierachy. (default: True)
        '''
        for Term_txt in obo_txt.split("[Term]\n"):
            if not Term_txt.strip():
                continue
            Term=GO_Term(Term_txt)
            namespace=Term.namespace
            if not Term.namespace in GO_namespace_to_Aspect:
                sys.stderr.write("ERROR! Unkown namespace %s for %s"%(
                    Term.namespace, Term.id))
                continue
            else:
                Aspect=GO_namespace_to_Aspect[namespace]

            if "is_obsolete: true" in Term_txt: # obsolete GO Term
                self[Aspect]["is_obsolete"].append(Term.id)
                if Term.id in self[Aspect]["Term"]:
                    del self[Aspect]["Term"][Term.id]
            else: # non-obsolete GO Term
                self[Aspect]["Term"][Term.id]=Term

            if uninformative_txt in Term_txt: # uninformative GO Term
                self[Aspect]["uninformative"].append(Term.id)

            for alt_id in Term.alt_id:
                self[Aspect]["alt_id"][alt_id]=Term.id
        if update_hierachy:
            self.update_is_a()
        return

    def update_is_a(self,Aspect=''):
        '''update self[Aspect]["is_a"] GO hierachy'''
        if not Aspect: # update all Aspect
            for Aspect in self:
                self.update_is_a(Aspect)
            return

        ## update direct parent ##
        parent_level=0
        add_parent=False
        print >>sys.stderr,"updating parent level %u for %s"%(parent_level,Aspect)
        for Term_id,Term in self[Aspect]["Term"].items():
            if Term.is_a:
                if not Term_id in self[Aspect]["is_a"]:
                    self[Aspect]["is_a"][Term_id]=set(Term.is_a)
                    add_parent=True
                elif Term.is_a - self[Aspect]["is_a"][Term_id]:
                    self[Aspect]["is_a"][Term_id] |=Term.is_a
                    add_parent=True
        ## indirect parent ##
        while(add_parent): # loop while new parent is added
            parent_level+=1
            add_parent=False
            print >>sys.stderr,"updating parent level %u for %s"%(parent_level,Aspect)
            for Term_id,parent_set in self[Aspect]["is_a"].items():
                indirect_parent_set=set()
                for parent in parent_set:
                    parent_id=parent[:10]
                    if parent_id in self[Aspect]["is_a"]:
                        indirect_parent_set|=self[Aspect]["is_a"][parent_id]

                if indirect_parent_set - self[Aspect]["is_a"][Term_id]:
                    add_parent=True
                    self[Aspect]["is_a"][Term_id] |= indirect_parent_set
        return

    def obsolete(self):
        '''return a string listing all obsolete Term'''
        return '\n'.join(['\n'.join(self[Aspect]["is_obsolete"]) \
            for Aspect in self])+'\n'
    
    def uninformative(self,name=True):
        '''return a string listing all uninformative Term
        name - True:  (default) show name alongside GO id
               False: only show GO id
        '''
        if name==False:
            uninformative_str='\n'.join(['\n'.join(self[Aspect]["uninformative"]) \
            for Aspect in self])+'\n'
        else:
            uninformative_str=''
            for Aspect in self:
                for Term_id in self[Aspect]["uninformative"]:
                    uninformative_str+=self.short(Term_id)
        return uninformative_str

    def Term(self,Term_id=""):
        '''return a GO_Term class for a specific Term
        or a list of Term belong to an Aspect'''
        if Term_id.startswith("GO:"):
            Term_id=self.alt_id(Term_id)
            for Aspect in self:
                if Term_id in self[Aspect]["Term"]:
                    return self[Aspect]["Term"][Term_id]
        elif not Term_id:
            return self.Term('F')+self.Term('P')+self.Term('C')
        elif Term_id in GO_namespace_to_Aspect and \
            GO_namespace_to_Aspect[Term_id] in self:
            Aspect=GO_namespace_to_Aspect[Term_id]
            return [v for k,v in self[Aspect]["Term"].items()]

        sys.stderr.write("ERROR! Cannot find GO Term %s\n"%Term_id)
        return []

    def has_a(self,Term_id=""):
        '''cumulatively list all child GO id'''
        sys.stderr.write("ERROR! Unimplemented.\n")
        return

    def alt_id(self,Term_id=""):
        '''return primary id of an alt_id
        If the input is primary id, return itself
        If the input is not found, return empty
        '''
        if Term_id.startswith("GO:"):
            for Aspect in self:
                if Term_id in self[Aspect]["Term"]:
                    return Term_id
                elif Term_id in self[Aspect]["alt_id"]:
                    return self[Aspect]["alt_id"][Term_id]
        elif not Term_id:
            return ''.join([self.alt_id(Aspect) for Aspect in self])
        elif Term_id in GO_namespace_to_Aspect and \
            GO_namespace_to_Aspect[Term_id] in self:
            Aspect=GO_namespace_to_Aspect[Term_id]
            return '\n'.join([k+'   '+v for k,v in \
                self[Aspect]["alt_id"].items()])+'\n'

        sys.stderr.write("ERROR! Cannot find GO Term '%s'"%Term_id)
        return ''

    def is_a(self,Term_id='',direct=True,name=False,number=False):
        '''List all parent GO id

        direct - True:  (default) only direct parent GO or
                 False: cumulatively list all parent GO
        name   - True:  (default) show name alongside GO id
                 False: only show GO id
        number - False: (default) list parent GO id
                 True: only list parent GO number
        '''
        Term_list=self.Term(Term_id)
        Term_list=[Term_list] if not isinstance(Term_list,list) else Term_list

        is_a_list=[]
        for Term in Term_list:
            is_a_entry=[Term.short() if name else Term.id]
            if direct:
                is_a_entry+=[is_a_id if name else is_a_id[:10]
                    for is_a_id in Term.is_a]
            else: # cumulatively list all parent
                for Aspect in self:
                    if Term.id in self[Aspect]["is_a"]:
                        is_a_entry+=[is_a_id if name else is_a_id[:10]
                        for is_a_id in self[Aspect]["is_a"][Term.id] ]
            is_a_list.append((len(is_a_entry)-1,is_a_entry))

        is_a_list=sorted(is_a_list)
        is_a_list=[is_a_entry[0]+'\t%u'%is_a_number if number else \
            '\t'.join(is_a_entry) for is_a_number,is_a_entry in is_a_list]
        return '\n'.join(is_a_list)+'\n'

    def csv(self,Term_id=''):
        '''return text for GO Term in CSV format'''
        Term_list=self.Term(Term_id)
        Term_list=[Term_list] if not isinstance(Term_list,list) else Term_list
        return '\n'.join([Term.__str__() for Term in Term_list])+'\n'

    def obo(self,Term_id=''):
        '''return text for GO Term in OBO format'''
        Term_list=self.Term(Term_id)
        Term_list=[Term_list] if not isinstance(Term_list,list) else Term_list
        return '\n'.join([Term.__repr__() for Term in Term_list])+'\n'
    
    def short(self,Term_id=''):
        '''return text for GO Term in short format: id ! name'''
        Term_list=self.Term(Term_id)
        Term_list=[Term_list] if not isinstance(Term_list,list) else Term_list
        return '\n'.join([Term.short() for Term in Term_list])+'\n'

    def summary(self,Aspect=''):
        '''for Aspect "Aspect", return a list summarizing number of
        Term, alt_id, is_obsolete'''
        Term_num=0
        is_obsolete_num=0
        alt_id_num=0
        if Aspect: # summary for one Aspect
            Term_num=len(
                self[GO_namespace_to_Aspect[Aspect]]["Term"])
            alt_id_num=len(
                self[GO_namespace_to_Aspect[Aspect]]["alt_id"])
            is_obsolete_num=len(
                self[GO_namespace_to_Aspect[Aspect]]["is_obsolete"])
        else:     # summary for all three Aspect
            for Aspect in self:
                Term_num_Aspect,alt_id_num_Aspect,is_obsolete_num_Aspect= \
                    self.summary(Aspect)
                Term_num       +=Term_num_Aspect
                is_obsolete_num+=is_obsolete_num_Aspect
                alt_id_num     +=alt_id_num_Aspect
        return [Term_num,alt_id_num,is_obsolete_num]
    
    def __str__(self):
        '''short summary for the whole obo class'''
        summary_txt='Aspect\tTerm\talt_id\tis_obsolete\n'
        for Aspect in self:
            summary_txt+='%s\t%s\n'%(Aspect,
                '\t'.join(map(str,self.summary(Aspect))))
        summary_txt+='%s\t%s\n'%("all",
                '\t'.join(map(str,self.summary())))
        return summary_txt

def parse_obo_txt(obo_txt=''):
    '''parse OBO format gene ontogology plain text "obo_txt", 
    return an obo class'''
    obo_txt=obo_txt[obo_txt.find("[Term]")-1:]
    obo_txt=obo_txt[:obo_txt.find("[Typedef]")]
    obo_dict=obo(obo_txt)
    return obo_dict

def parse_obo_file(obo_file="go-basic.obo"):
    '''parse OBO format gene ontogology file "obo_file", 
    return an obo class'''
    fp=open(obo_file,'rU')
    obo_txt=fp.read()
    fp.close()
    obo_dict=parse_obo_txt(obo_txt)
    return obo_dict

def obo2csv(obo_file="go-basic.obo",prefix=''):
    '''convert obo format gene ontogology file "obo_file" into CSV format 
    tabular files (all files prefixed by "prefix"):
    go.obo.F.csv,   go.obo.P.csv,   go.obo.C.csv
    go.obo.F.list,  go.obo.P.list,  go.obo.C.list
    go.obo.F.is_a,  go.obo.P.is_a,  go.obo.C.is_a
    go.obo.F.has_a, go.obo.P.has_a, go.obo.C.has_a
    go.obo.obsolete
    go.obo.GOaltids.dat
    go.obo.uninformative
    '''
    obo_dict=parse_obo_file(obo_file)
    file_list=[]
    basename=os.path.basename(obo_file)

    for Aspect in obo_dict:
        # *.csv contains
        # id, name, Aspect, is_a (direct parent)
        filename=prefix+basename+'.'+Aspect+".csv"
        fp=open(filename,'w')
        fp.write(obo_dict.csv(Aspect))
        fp.close()
        file_list.append(os.path.abspath(filename))
        
        # *.list contains
        # id ! name
        filename=prefix+basename+'.'+Aspect+".list"
        fp=open(filename,'w')
        fp.write(obo_dict.short(Aspect))
        fp.close()
        file_list.append(os.path.abspath(filename))

        # *.is_a contains
        # id\tis_a (culmulative)
        filename=prefix+basename+'.'+Aspect+".is_a"
        fp=open(filename,'w')
        fp.write(obo_dict.is_a(Aspect,direct=False,name=True))
        fp.close()
        file_list.append(os.path.abspath(filename))

    # *.obsolete lists obsolete Terms
    filename=prefix+basename+".obsolete"
    fp=open(filename,'w')
    fp.write(obo_dict.obsolete())
    fp.close()
    file_list.append(os.path.abspath(filename))

    # *.uninformative lists obsolete Terms
    filename=prefix+basename+".uninformative"
    fp=open(filename,'w')
    fp.write(obo_dict.uninformative(name=True))
    fp.close()
    file_list.append(os.path.abspath(filename))

    # *.GOaltids.dat contains:
    #alt_id   id
    filename=prefix+basename+".GOaltids.dat"
    fp=open(filename,'w')
    fp.write(obo_dict.alt_id())
    fp.close()
    file_list.append(os.path.abspath(filename))

    print obo_dict.__str__()
    return file_list

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    file_list=obo2csv(sys.argv[1])
    print "output files:"
    for filename in file_list:
        print filename
