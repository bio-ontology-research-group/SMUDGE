@Grapes([
	  @Grab(group='org.semanticweb.elk', module='elk-owlapi', version='0.4.3'),
	  @Grab(group='net.sourceforge.owlapi', module='owlapi-api', version='4.2.5'),
	  @Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='4.2.5'),
	  @Grab(group='net.sourceforge.owlapi', module='owlapi-impl', version='4.2.5'),
	  @Grab(group='net.sourceforge.owlapi', module='owlapi-parsers', version='4.2.5'),
	  @Grab(group='org.apache.jena', module='apache-jena-libs', version='3.1.0', type='pom')
	])

import org.semanticweb.owlapi.model.parameters.*
import org.semanticweb.elk.owlapi.ElkReasonerFactory;
import org.semanticweb.elk.owlapi.ElkReasonerConfiguration
import org.semanticweb.elk.reasoner.config.*
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.reasoner.*
import org.semanticweb.owlapi.reasoner.structural.StructuralReasoner
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary;
import org.semanticweb.owlapi.model.*;
import org.semanticweb.owlapi.io.*;
import org.semanticweb.owlapi.owllink.*;
import org.semanticweb.owlapi.util.*;
import org.semanticweb.owlapi.search.*;
import org.semanticweb.owlapi.manchestersyntax.renderer.*;
import org.semanticweb.owlapi.reasoner.structural.*
import org.apache.jena.rdf.model.*
import org.apache.jena.util.*

def cli = new CliBuilder()
cli.with {
usage: 'Self'
  h longOpt:'help', 'this information'
  i longOpt:'input', 'input RDF file', args:1, required:true
  u longOpt:'undirected', 'build undirected graph (default: false)', args:1, required:false
  c longOpt:'classify', 'use an OWL reasoner to classify the RDF dataset (must be in RDF/XML) before graph generation (default: false)', args:1, required:false
  f longOpt:'format', 'RDF format; values are "RDF/XML", "N-TRIPLE", "TURTLE" and "N3" (default: RDF/XML)', args:1, required:false
  d longOpt:'ontology-directory', 'directory with ontologies to use for reasoning', args:1, required:false
  o longOpt:'output', 'output corpus file',args:1, required:true

}
def opt = cli.parse(args)
if( !opt ) {
  //  cli.usage()
  return
}
if( opt.h ) {
    cli.usage()
    return
}

def undirected = false
if (opt.u && opt.u != "false") {
  undirected = true
}
def classify = false
if (opt.c && opt.c != "false") {
  classify = true
}
def format = "RDF/XML"
if (opt.f) {
  format = opt.f
}



def f = File.createTempFile("temp",".tmp")
if (classify) {
  OWLOntologyManager manager = OWLManager.createOWLOntologyManager()
  def oset = new LinkedHashSet()
  oset.add(manager.loadOntologyFromOntologyDocument(new File(opt.i)))
  if (opt.d) {
    new File(opt.d).eachFile { ofile ->
      oset.add(manager.loadOntologyFromOntologyDocument(ofile))
    }
  }
  OWLOntology ont = manager.createOntology(IRI.create("http://aber-owl.net/rdfwalker/t.owl"),oset)
  OWLDataFactory fac = manager.getOWLDataFactory()
  ConsoleProgressMonitor progressMonitor = new ConsoleProgressMonitor()
  OWLReasonerConfiguration config = new SimpleConfiguration(progressMonitor)
  ElkReasonerFactory f1 = new ElkReasonerFactory()
  OWLReasoner reasoner = f1.createReasoner(ont,config)
  def cc = 0
  new InferredClassAssertionAxiomGenerator().createAxioms(fac, reasoner).each { ax ->
    manager.addAxiom(ont, ax)
    cc += 1
  }
  println "$cc axioms inferred."

  PrintWriter fout = new PrintWriter(new BufferedWriter(new FileWriter(opt.o)))

  //walk ontology upward dfs for each data node
  new File("../../Documents/smudge_data/pheno2_data.txt").splitEachLine("\t") {line ->
  def node = line[0]
  def class_iri = line[1]
  def curr_path = [node,class_iri]
  def classes_stack = []
  classes_stack.push(curr_path)

  if (ont.containsClassInSignature(IRI.create(class_iri))){
   while(classes_stack){
     curr_path = classes_stack.pop()
     last_class = curr_path[-1]

   if (last_class.indexOf('Thing') > -1){
      fout.print(curr_path[0] + '\t' + curr_path[1].split('/')[-1] + '\t')
     for (item in curr_path[2..-1]){
      shortitem = item.split('/')[-1] 
      if(shortitem.startsWith('HP_') || shortitem.startsWith('MP_')){ //select only HP and MP classes from PhenomNet
          fout.print('subClassOf'+'\t'+ shortitem +'\t')
        }
      }
      
      fout.print('\n')
      continue
    }
   OWLClass c = fac.getOWLClass(IRI.create(last_class))
   def superClasses = reasoner.getSuperClasses(c, true)
   for (Node<OWLClass> parent: superClasses){
     def superclass = parent.getRepresentativeElement().toString().replaceAll(">","").replaceAll("<","")
     new_path = curr_path + [superclass]
     classes_stack.push(new_path)
    }   
   }
  }
 }
 fout.flush()
 fout.close()
}