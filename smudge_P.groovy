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
    
  manager.saveOntology(ont, IRI.create(f.toURI()))
  println "$cc axioms inferred." 


PrintWriter fout = new PrintWriter(new BufferedWriter(new FileWriter(opt.o)))

new File("../../Documents/smudge_data/pheno2_data.txt").splitEachLine("\t") {line ->
def node = line[0]
def class_iri = line[1]
fout.print(node+ ' ')
if (ont.containsClassInSignature(IRI.create(class_iri))){
  OWLClass c = fac.getOWLClass(IRI.create(class_iri))
  NodeSet superclasses = reasoner.getSuperClasses(c, false);
  fout.print(c.getIRI().getShortForm() + " ")
  for (Node<OWLClass> superc: superclasses){
    fout.print(superc.getRepresentativeElement().getIRI().getShortForm() + " ")
      }
  }
  fout.print('\n')
}

fout.flush()
fout.close()
}