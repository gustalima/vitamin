#!/usr/bin/env python
# -*- coding: utf-8 -*-

from gi.repository import Gtk
from gi.repository import Gdk
from gi.repository import GObject
from gi.repository import GLib
from collections import Counter
import pymol
import threading
import thread
import shutil
import os
import glob
import time
import math
import modeller
from time import sleep
from threading import Timer
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW
from Bio.PDB import *
from gi.repository.Pango import FontDescription
from Bio.Blast import NCBIXML
from multiprocessing import Process


GObject.threads_init()
GLib.threads_init()
Gdk.threads_init()


initfiles = Counter(glob.glob('*.*'))


class FuncThread(threading.Thread):
        def __init__(self, target, *args):
                self._target = target
                self._args = args
                self._stop = threading.Event()
                threading.Thread.__init__(self)

        def run(self):
                self._target(*self._args)

        def stop(self):
                self._stop.set()

class Handler():
        def __init__(self):
                global parameters
                global milestones
                # try:
                #     conf = open('config','r').read().split('\n')
                #     parameters = conf[0].split(' ')
                #     milestones = conf[1].split(' ')
                # except:
                #     parameters = ['fasta','1pdb','2pdb','3pdb','1chain','2chain','3chain', 'alifile']
                #     milestones = ['align','readseq','model', 'bestmodel', 'ramachandran', 'evaluation', 'ref_bestmodel']

                parameters = ['fasta','1pdb','2pdb','3pdb','1chain','2chain','3chain', 'alifile']
                milestones = ['align','readseq','model', 'bestmodel', 'ramachandran', 'evaluation', 'ref_bestmodel']

                builder.get_object('button10').set_sensitive(True)
                builder.get_object('button49').set_sensitive(True)
                self.create_treeview()

#       XML file filter for dialogs
                self.xmlfilter = Gtk.FileFilter()
                self.xmlfilter.add_pattern('*.xml')

#       Combobox for blast algorithm selection
                self.blastengine = builder.get_object('comboboxtext1')
                self.cell = Gtk.CellRendererCombo()
                self.blastalgo = ['Blastp', 'PSI-BLAST', 'PHI-BLAST', 'DELTA-BLAST']
                for each in self.blastalgo:
                        self.blastengine.append_text(each)
                self.blastengine.set_active(0)

#       Textview for protein sequence
                self.textview = builder.get_object('textview2')
                self.textview.modify_font(FontDescription("sans 8"))



        def on_button47_clicked(self,*args):
                '''Checks and calls a function for Load XML Blast output from previous search'''
                self.control = self.get_proteinsequence()
                if self.control != "OK":
                        return
                else:
                        pass
                self.choosefile = Gtk.FileChooserDialog(title = "Select the Blast XML outputfile",
                                                        action = Gtk.FileChooserAction.OPEN,
                                                        buttons = ["Choose", Gtk.ResponseType.OK,
                                                                   "Cancel", Gtk.ResponseType.CANCEL])

                self.choosefile.set_filter(self.xmlfilter)
                self.response = self.choosefile.run()
                if self.response == Gtk.ResponseType.OK:
                        self.blastresultfile = self.choosefile.get_filename()
                        self.choosefile.destroy()
                        self.load_blastresults(self.blastresultfile)

                elif self.response == Gtk.ResponseType.CANCEL:
                        pass
                self.choosefile.destroy()

        def on_textview2_key_release_event(self,*args):
                '''Makes the protein sequence clean and tight'''
                self.buffer = self.textview.get_buffer()
                self.startiter = self.buffer.get_start_iter()
                self.enditer = self.buffer.get_end_iter()
                self.text_buffer = str(self.buffer.get_text(self.startiter, self.enditer, True))
                self.text = self.text_buffer.replace('\n', '')
                self.text = self.text.replace(' ', '')
                self.buffer.set_text(self.text)

        def get_proteinsequence(self, *args):
                '''Performs some checks of consistency of protein sequence and takes the sequence
                from the Textview'''
                self.on_textview2_key_release_event()
                self.buffer = self.textview.get_buffer()
                self.startiter = self.buffer.get_start_iter()
                self.enditer = self.buffer.get_end_iter()
                self.proteinsequence = str(self.buffer.get_text(self.startiter, self.enditer, True))
                self.sequencelength = len(self.proteinsequence)
                if '>' in self.proteinsequence:
                        builder.get_object('label41').set_text('You must enter just the aminoacid sequence \n without any special character')
                        builder.get_object('warning_dialog1').show()
                        self.buffer.set_text('')
                        return
                elif self.sequencelength == 0:
                        builder.get_object('label41').set_text('You must enter the aminoacid sequence')
                        builder.get_object('warning_dialog1').show()
                        self.buffer.set_text('')
                        return
                elif any(c.isdigit() for c in self.proteinsequence):
                        builder.get_object('label41').set_text('You must enter a valid aminoacid sequence')
                        builder.get_object('warning_dialog1').show()
                        self.buffer.set_text('')
                        return
                else:
                        return "OK"

        def on_button41_clicked(self, *args):
                '''Runs the Blast Search with the inputed protein sequence'''
                self.control = self.get_proteinsequence()
                if self.control != "OK":
                        return
                else:
                        self.start_spinner("Blasting sequence")

                def thread_run2(*args):
                        '''Thread for all blasting process'''
                        sleep(0.5)
                        run_control = blast_sequence()
                        if run_control == 'timeout':
                                Gdk.threads_enter()
                                self.spinnerdialog.destroy()
                                builder.get_object('label41').set_text('Timeout! \n There was a problem during blast job. \n Try again latter! ')
                                builder.get_object('warning_dialog1').show()
                                Gdk.threads_leave()
                        else:
                                Gdk.threads_enter()
                                self.spinnerdialog.destroy()
                                self.load_blastresults('my_blast.xml')
                                Gdk.threads_leave()
                        return

                def blast_helper(proteinseq):
                        '''Helper function that perform the blast job '''
                        blastresults = NCBIWWW.qblast('blastp', 'pdb', proteinseq)
                        save_file = open("my_blast.xml", 'w')
                        save_file.write(blastresults.read())

                def blast_sequence():
                        '''Sets blasting job at another process and checks for the timeout'''
                        p = Process(target=blast_helper, args=(self.proteinsequence,))
                        p.start()
                        self.start_time = time.time()
                        while time.time() - self.start_time < 300:
                                if p.is_alive():
                                        sleep(10)
                                        pass
                                else:
                                        p.terminate()
                                        return
                        p.terminate()
                        return 'timeout'

                t = threading.Thread(target=thread_run2)
                t.start()

        def load_blastresults(self, blastfile):
                '''Loads results from Blast Search '''
                self.blastfile = open(blastfile, 'r')
                self.blastrecords = NCBIXML.read(self.blastfile)
                self.treeview = builder.get_object('treeview1').get_model().clear()

                for x, y in zip(self.blastrecords.alignments, self.blastrecords.descriptions):
                        for each in x.hsps:
                                self.coverage = (each.query_end - each.query_start)
                                self.coveragepercent = str(round(self.coverage / float(self.sequencelength)*100))+' %'
                                self.identity = str((round(each.identities / float(self.coverage)*100)))+' %'
                                self.evalue = str("%.1g" % each.expect)
                        self.pdb_id = y.title.split('|')[3]
                        self.template_descripton = y.title.split('>')[0]
                        if float(self.identity.split(' %')[0]) >= 40.0:
                                self.backgroundcolor = '#00FF00'
                        elif 25.0 <= float(self.identity.split(' %')[0]) < 40.0:
                                self.backgroundcolor = '#FFFF00'
                        elif 15.0 <= float(self.identity.split(' %')[0]) < 25.0:
                                self.backgroundcolor = '#FF9933'
                        elif 0.0 <= float(self.identity.split(' %')[0]) < 15.0:
                                self.backgroundcolor = '#FF0000'


                        self.treeview_liststore.append([self.backgroundcolor,
                                                        self.template_descripton,
                                                        self.pdb_id,
                                                        self.identity,
                                                        self.coveragepercent,
                                                        self.evalue])


        def start_spinner(self,label):
                '''Starts spinner dialog with a flag for waiting for the job (label) to be Done!'''
                self.spinnerdialog = Gtk.Dialog(title = "Please wait!")
                self.spinnerdialog.set_default_size(300,150)
                self.spinnerdialog.set_decorated(False)
                self.spinnerlabel = Gtk.Label(label)
                self.spinnercontent = self.spinnerdialog.get_content_area()
                self.spinner = Gtk.Spinner()
                self.spinner.show()
                self.spinnercontent.pack_end(self.spinnerlabel, False, False, 10)
                self.spinnercontent.pack_start(self.spinner, True, True, 10)
                self.spinner.start()
                self.spinnerlabel.show()
                self.spinnerdialog.set_position(Gtk.WindowPosition.CENTER)
                self.spinnerdialog.show()
                return


        def create_treeview(self):
                '''Creates the treeview for blast results output'''
                self.treeview = builder.get_object('treeview1')
                self.treeview_liststore = Gtk.ListStore(str, str, str, str, str, str)
                self.treeview.set_model(self.treeview_liststore)
                self.treeview.set_grid_lines(True)
                self.cell_rend = Gtk.CellRendererText()
                self.cell_rend2 = Gtk.CellRendererText()
                self.pango = FontDescription("sans 8")
                self.cell_rend2.set_property('font-desc', self.pango)
                self.cell_rend.set_alignment(0.5,0.5)
                self.cell_rend2.set_alignment(0.0,0.5)
                self.cell_rend2.props.wrap_width = 600
                self.cell_rend2.props.wrap_mode = Gtk.WrapMode.WORD
                self.selection_source = self.treeview.get_selection()
                self.selection_source.set_mode(Gtk.SelectionMode.MULTIPLE)

                # SETTING COLUMS TO THE TREEVIEW
                ######STATUS COLUMN######
                self.column1 = Gtk.TreeViewColumn("S",Gtk.CellRendererText(), background=0, )
                self.column1.set_resizable(False)
                self.column1.set_alignment(0.5)
                self.column1.set_expand(False)
                self.column1.set_clickable(True)


                ######DESCRIPTION COLUMN######
                self.column2 = Gtk.TreeViewColumn("Description",self.cell_rend2, text=1)
                self.column2.set_resizable(False)
                self.column2.set_alignment(0.5)
                self.column2.set_expand(True)
                self.column2.set_max_width(50)
                self.column2.set_sort_column_id(1)

                ######PDB ID COLUMN######
                self.column3 = Gtk.TreeViewColumn("PDB ID",self.cell_rend, text=2)
                self.column3.set_resizable(False)
                self.column3.set_alignment(0.5)
                self.column3.set_expand(False)
                self.column3.set_sort_column_id(2)

                ######IDENTITY COLUMN######
                self.column4 = Gtk.TreeViewColumn("Identity",self.cell_rend, text=3)
                self.column4.set_resizable(False)
                self.column4.set_alignment(0.5)
                self.column4.set_expand(False)
                self.column4.set_sort_column_id(3)

                ######Coverage COLUMN######
                self.column5 = Gtk.TreeViewColumn("Coverage",self.cell_rend, text=4)
                self.column5.set_resizable(False)
                self.column5.set_alignment(0.5)
                self.column5.set_expand(False)
                self.column5.set_sort_column_id(4)

                ######E VALUE COLUMN######
                self.column6 = Gtk.TreeViewColumn("e-val",self.cell_rend, text=5)
                self.column6.set_resizable(False)
                self.column6.set_alignment(0.5)
                self.column6.set_expand(False)


                self.treeview.insert_column(self.column1, -1)
                self.treeview.insert_column(self.column2, -1)
                self.treeview.insert_column(self.column3, -1)
                self.treeview.insert_column(self.column4, -1)
                self.treeview.insert_column(self.column5, -1)
                self.treeview.insert_column(self.column6, -1)

        def on_button10_clicked(self,*args):
                self.treeview = builder.get_object('treeview1')
                ''' Fetchs the selected PDB from blast search output'''
                self.model, self.selected_rows = self.treeview.get_selection().get_selected_rows()
                if self.selected_rows == []:
                        return
                else:
                        self.pdbcodes = []
                        for each in self.selected_rows:
                                self.pdbcodes.append(self.model[each][2])
                        pymol.cmd.do('delete all; fetch %s, async=0; split_chains; as cartoon; util.cbc; alignto %s_A; zoom %s_A' %(' '.join(self.pdbcodes),self.pdbcodes[0],self.pdbcodes[0]))
                        while len(self.pdbcodes)!= 3:
                                self.pdbcodes.append('')
                        parameters[1:4] = self.pdbcodes


        def on_button49_clicked(self, *args):
                '''General function for retrieving additional information of pdb entries'''
                self.model, self.selected_rows = self.treeview.get_selection().get_selected_rows()
                if len(self.selected_rows) == []:
                        return
                else:
                        self.start_spinner("Retrieving pdb information")

                t = threading.Thread(target=self.retrieve_information)
                t.start()

        def retrieve_information(self,*args):
                sleep(0.5)
                run_control = self.retrieve_helper()
                if run_control == 'timeout':
                        Gdk.threads_enter()
                        self.spinnerdialog.destroy()
                        builder.get_object('label41').set_text('Timeout! \n There was a problem trying to get pdb information. \n Try again latter! ')
                        builder.get_object('warning_dialog1').show()
                        Gdk.threads_leave()
                else:
                        Gdk.threads_enter()

                        for each in self.selected_rows:
                                self.pdbadress = './pdbs/'+str(self.model[each][2]) + '.pdb'
                                self.pdbfile = open(self.pdbadress,'r')
                                self.journalinf = ''
                                self.refinement_r = ''
                                self.completeness = ''
                                self.refinement_rfree = ''
                                self.completeness_helper = 0
                                self.rvalue_helper = 0
                                self.rfree_helper = 0
                                for each2 in self.pdbfile:

                                        self.keyvalue = each2.split()[0]
                                        if self.keyvalue == "JRNL":
                                                if each2.split()[1] == 'REF':
                                                        for x in each2.split()[2:]:
                                                                self.journalinf = self.journalinf + ' ' + str(x)
                                        elif self.keyvalue == 'REMARK':
                                                if "COMPLETENESS" in each2:
                                                        if self.completeness_helper == 0:
                                                                self.completeness = self.completeness + ' ' + str(each2.split()[-1])
                                                                self.completeness_helper += 1
                                                                self.completeness
                                                        else:
                                                                pass
                                                if 'R VALUE' and 'WORKING SET' in each2:
                                                        if self.rvalue_helper == 0:
                                                                self.refinement_r = self.refinement_r + ' ' + str(each2.split()[-1])
                                                                self.rvalue_helper+= 1
                                                                print self.refinement_r
                                                        else:
                                                                pass
                                                if 'FREE R VALUE' in each2:
                                                        if "TEST" in each2:
                                                                pass
                                                        else:
                                                                if self.rfree_helper == 0:
                                                                        self.refinement_rfree = self.refinement_rfree + ' ' + str(each2.split()[-1])
                                                                        self.rfree_helper += 1
                                                                        print self.refinement_rfree
                                                                else:
                                                                        pass
                                self.pdbfile.close()
                                self.pdbfile = open(self.pdbadress,'r')
                                self.header_dict=parse_pdb_header(self.pdbfile)
                                print self.header_dict['name']
                                print self.header_dict['source']
                                self.model[each][1] = ("Title: " + str(self.header_dict['name']) + '\n' +
                                                        "Structure method: " + str(self.header_dict['structure_method']) + '\n' +
                                                        "Source: " + str(self.header_dict['source']['1']['organism_scientific']) + '\n' +
                                                        "Resolution: " + str(self.header_dict['resolution']) + '\n' +
                                                        "Deposition date: " + str(self.header_dict['deposition_date']) + '\n' +
                                                        "Authors: " + str(self.header_dict['author']) + '\n' +
                                                        'Journal ref:' + self.journalinf + '\n\n' +
                                                        "Refinement details:" + '\n'+
                                                        "Completeness:" + self.completeness + '\n' +
                                                        'Rvalue:' + self.refinement_r + '\n' +
                                                        'FreeRvalue:' + self.refinement_rfree+ '\n\n')
                                self.pdbfile.close()

                        self.pdbfile.close()
                        self.spinnerdialog.destroy()
                        Gdk.threads_leave()
                return




        def retrieve_helper(self,*args):
                '''Calls the retrieve pdb info funcion helper using another process'''
                p = Process(target=self.retrieve_worker)
                p.start()
                self.start_time = time.time()
                while time.time() - self.start_time < 300:
                        if p.is_alive():
                                sleep(5)
                                pass
                        else:
                                p.terminate()
                                return
                p.terminate()
                return 'timeout'

        def retrieve_worker(self,*args):
                '''Downloads selected pdb files and extracts additional information'''
                self.parser = PDBParser(PERMISSIVE=1)
                self.pdblist = PDBList()
                for each in self.selected_rows:
                        self.pdblist.retrieve_pdb_file(self.model[each][2], pdir='./pdbs')
                        self.oldname = glob.glob('./pdbs/*.ent')[0]
                        self.newname = "./pdbs/"+self.model[each][2]+'.pdb'
                        os.system('mv %s %s' %(self.oldname, self.newname))
                return


        def onDeleteWindow(self, *args):
                Gtk.main_quit(*args)

        def on_main_destroy(self, *args):
                open('config','w').write(','.join(parameters)+'\n'+','.join(milestones))
                Gtk.main_quit(*args)

        def exit(widget, callback_data=None):           #quit for menubar
                builder.get_object('main').destroy()
                Gtk.main_quit()

        def save_workspace(*args):


                endfiles = Counter(glob.glob('*.*'))
                savefiles = list(endfiles-initfiles)
                directory = "Results_"
                directory += time.ctime().replace(' ', '_')


                try:                                    #create folder for results
                        os.stat(directory)
                except:
                        os.mkdir(directory)

                for files in savefiles:                 #move every file create since last session to results
                        shutil.move(files, directory)
                if not os.stat('Results/'):
                        os.mkdir('Results')
                for folder in glob.glob('Results_*'):
                                shutil.move(folder, 'Results/')


#       Window open/close instructions for PDB files

        def open_dialog1(self, *args):
                builder.get_object('filechooserdialog1').run()
                builder.get_object('filechooserdialog1').hide()

        def open_dialog2(self, *args):
                builder.get_object('filechooserdialog2').run()
                builder.get_object('filechooserdialog2').hide()

        def open_working_directory(self, *args):
                builder.get_object('filechooserdialog3').run()
                builder.get_object('filechooserdialog3').set_action(gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER)
                builder.get_object('filechooserdialog3').hide()

        def cancel_working_directory(self, *args):
                builder.get_object('filechooserdialog3').hide()

        def choose_working_directory(self, *args):
                working_dir = builder.get_object('filechooserdialog3').get_filenames()
                builder.get_object('filechooserdialog3').hide()
                os.system('lsyncd -log all -delay 2 -nodaemon -rsync %s %s' %(os.getcwd(), working_dir[0]))

        def reinit_vitamin(self, *args):
                pymol.cmd.do('delete all')
                os.system('python reinit.py')

        def cancel_dialog1(self, *args):
                builder.get_object('filechooserdialog1').hide()

        def cancel_dialog2(self, *args):
                builder.get_object('filechooserdialog2').hide()

        def close_warning_dialog1(self, *args):
                builder.get_object('warning_dialog1').hide()


        def get_profile(self, profile_file, seq):
            """Read `profile_file` into a Python array, and add gaps corresponding to
               the alignment sequence `seq`."""
            # Read all non-comment and non-blank lines from the file:
            f = file(profile_file)
            vals = []
            for line in f:
                if not line.startswith('#') and len(line) > 10:
                    spl = line.split()
                    vals.append(float(spl[-1]))
            # Insert gaps into the profile corresponding to those in seq:
            for n, res in enumerate(seq.residues):
                for gap in range(res.get_leading_gaps()):
                    vals.insert(n, None)
            # Add a gap at position '0', so that we effectively count from 1:
            vals.insert(0, None)
            return vals

        def create_alignment(self, *args):

                chaincodes = pymol.cmd.get_names(enabled_only=1)

                for pdb in chaincodes:
                        pymol.cmd.do("save %s.pdb, %s" %(pdb,pdb))

                for pdb in chaincodes:
                        os.system('python src/readseq.py %s' % pdb)
                print 'Seg files for selected PDB chain(s) created sucessfully'


                buffering = builder.get_object('textview2').get_buffer()
                startiter = buffering.get_start_iter()
                enditer = buffering.get_end_iter()
                proteinsequence = str(buffering.get_text(startiter, enditer, True))
                #print proteinsequence


                if len(proteinsequence)>5:

                        #fastaseq = '>Search sequence\n'+proteinsequence

                        with open("input", "w") as text_file:
                                text_file.write('>Search sequence\n'+proteinsequence)
                        with open('input.seg', 'w') as segfile:
                                segfile.write('>P1;target\nsequence:target::::::::\n'+proteinsequence+'*')


                        fasta = 'input'



                if len(chaincodes) == 1:
                        os.system('python src/align.py %s %s | tee align.log' % (chaincodes[0], fasta))
                        parameters[7] = 'input_alignment.ali'
                        milestones[0] = 'OK'
                        print 'Single template alignment created sucessfully!'

                elif len(chaincodes)>1:
                        mcodes = ''
                        for code in chaincodes:
                                mcodes += code+'/'
                        if mcodes[-1] == '/':
                                mcodes = mcodes[:-1]
                        parameters[7] = 'input_alignment.ali'
                        os.system('python src/salign.py %s | tee salign.log' % mcodes)
                        milestones[0] = 'OK'
                        #os.system('python src/align2d_mult.py | tee align2d_mult.log' )
                        print len(chaincodes),' templates used'
                        print 'Multiple template alignment created sucessfully!'
                milestones[1] = 'OK'

                while len(chaincodes) < 3:
                        chaincodes.append('')
                st = ''
                if len(chaincodes)>3:
                        for i in chaincodes[3:]:
                                st+= str(i)+'/'

                chaincode = [chaincodes[0],chaincodes[1],st]
                parameters[4:7] = chaincode
                 #aqui vai precisar mudar o jeito que os pdbs sao armazenados na
                 #lista milestones... Na verdade precisa ver isso no codigo todo
                 #que ta foda



        def cleanup_pymol(self,*args):
                pdbs = pymol.cmd.get_names('all')
                chaincodes = pymol.cmd.get_names(enabled_only=1)

                for pdb in chaincodes:
                        # pdbs.remove(pdb)
                        if pdb not in chaincodes:
                            pymol.cmd.do('delete %s' % pdb)
                if len(chaincodes) == 1:
                        chaincodes.append('')
                        chaincodes.append('')
                elif len(chaincodes) == 2:
                        chaincodes.append('')
                st = ''
                if len(chaincodes)>3:
                        for i in chaincodes[3:]:
                                st+= str(i)+'/'

                chaincode = [chaincodes[0],chaincodes[1],st]
                parameters[4:7] = chaincode


#       Toggling buttons for custom activity

        def toggle_button_own_script(self, *args):
                if builder.get_object('checkbutton3').get_active() and builder.get_object('radiobutton1').get_active():
                        builder.get_object('filechooserbutton1').set_sensitive(True)
                else:
                        builder.get_object('filechooserbutton1').set_sensitive(False)

        def toggle_button_own_script2(self, *args):
                if builder.get_object('checkbutton3').get_active() and builder.get_object('radiobutton2').get_active():
                        builder.get_object('button46').set_sensitive(True)
                else:
                        builder.get_object('button46').set_sensitive(False)

        def toggle_script(self, *args):
                if builder.get_object('checkbutton3').get_active():
                        if builder.get_object('radiobutton1').get_active():
                                builder.get_object('filechooserbutton1').set_sensitive(True)
                        elif builder.get_object('radiobutton2').get_active():
                                builder.get_object('button46').set_sensitive(True)
                else:
                        builder.get_object('button46').set_sensitive(False)
                        builder.get_object('filechooserbutton1').set_sensitive(False)

#       Process FASTA file into seg format
        #
        # def basic_FASTA_input(self, *args):
        #         fasta = builder.get_object('basic_fasta_file').get_filename()
        #         parameters[0] = builder.get_object('basic_fasta_file').get_filename().split("/")[-1].split('.')[0]




#       Control About dialog

        def open_about(self, *args):
                builder.get_object('aboutdialog1').run()
                builder.get_object('aboutdialog1').hide()


#       Input Editor
        def edit_alignments(self, *args):
                chaincodes = pymol.cmd.get_names(enabled_only=1)
                file_list = ['input_alignment.ali', 'input_alignment.pap', 'input.seg', 'input']
                for code in chaincodes:
                    file_list.append(code+'.seg')

                for inpfile in file_list:
                    if os.path.isfile(inpfile) == False:
                        file_list.remove(inpfile)

                if len(file_list):
                    os.system('gedit %s &>null' %' '.join(file_list))
                else:
                    print "You don't have any input file in this folder. Try to create a new alignment with you PyMOL selection"



        def edit_basic_script(self, *args):
                os.system('gedit src/mod_model.py' )
                print 'Basic model script opened in your default document viewer'

#       Modeling functions

        def model_summary(self, *args):
                chaincodes = pymol.cmd.get_names(enabled_only=1)
                pdbs = ''
                for code in chaincodes:
                        pdbs +=code+'\n'
                pdbs = pdbs[:-1]
                try:
                        len(builder.get_object('basic_fasta_file').get_filename())
                        fasta = 'From file'
                except:
                        fasta = 'From BLAST search input'

                builder.get_object('summary').show()

                builder.get_object('label28').set_text(pdbs)
                builder.get_object('label24').set_text(builder.get_object('entry7').get_text())

        def model_summary_cancel(self,*args):
                builder.get_object('summary').hide()

        def basic_model(self, *args):
                builder.get_object('summary').hide()
                chaincodes = pymol.cmd.get_names(enabled_only=1)
                number = builder.get_object('entry7').get_text()
                hetatm = builder.get_object('checkbutton4').get_active()
                nproc  = builder.get_object('entry4').get_text()
                if nproc == '':
                        nproc = 1
                import multiprocessing
                if int(nproc) > multiprocessing.cpu_count():
                    nproc = multiprocessing.cpu_count()
                if len(chaincodes) == 1:

                        self.start_spinner('Creating models')

                        def thread_run(*args):
                                if builder.get_object('checkbutton3').get_active():
                                        if builder.get_object('radiobutton1').get_active():
                                                script = builder.get_object('basic_fasta_file').get_filename()
                                                model_ret = os.system('python %s %s %s | tee model.log' %(script, number, chaincodes))
                                                print 'Using your own script from %s' %script
                                        elif builder.get_object('radiobutton2').get_active():
                                                model_ret = os.system('python src/mod_model.py %s %s | tee model.log' %(number, chaincodes))
                                else:
                                        if hetatm:
                                                model_ret = os.system('python src/model.py %s %s %s true| tee model.log' %(number, chaincodes,nproc))
                                        else:
                                                model_ret = os.system('python src/model.py %s %s %s false| tee model.log' %(number, chaincodes,nproc))


                                GObject.idle_add(cleanup_func, model_ret)

                        def cleanup_func(model_ret):
                                Gdk.threads_enter()
                                self.spinnerdialog.destroy()
                                Gdk.threads_leave()

                                self.after_model()

                        t = threading.Thread(target=thread_run)
                        t.start()
                        milestones[2] = 'OK'

                elif len(chaincodes) > 1:
                        chaincodes = pymol.cmd.get_names(enabled_only=1)
                        pdbs = ''
                        for code in chaincodes:
                                pdbs += code+'/'
                        if pdbs[-1]=='/':
                                pdbs = pdbs[:-1]
                        fout=open('log.txt','w')
                        fout.write(pdbs)
                        fout.close()

                        self.start_spinner('Creating models')

                        def thread_run(*args):
                                if builder.get_object('checkbutton3').get_active():
                                        if builder.get_object('radiobutton1').get_active():
                                                script = builder.get_object('basic_fasta_file').get_filename()
                                                model_ret = os.system('python %s %s %s | tee model.log' %(script, number, pdbs))
                                                print 'Using your own script from %s' %script
                                        elif builder.get_object('radiobutton2').get_active():
                                                model_ret = os.system('python src/mod_multi_model.py %s %s | tee model.log' %(number, pdbs))

                                else:
                                        if hetatm:
                                                model_ret = os.system('python src/multi_model.py %s %s %s true| tee model.log' %(number, pdbs,nproc))
                                        else:
                                                model_ret = os.system('python src/multi_model.py %s %s %s false| tee model.log' %(number, pdbs,nproc))



                                GObject.idle_add(cleanup_func, model_ret)

                        def cleanup_func(model_ret):
                                Gdk.threads_enter()
                                self.spinnerdialog.destroy()
                                Gdk.threads_leave()

                                self.after_model()

                        t = threading.Thread(target=thread_run)
                        t.start()
                        milestones[2] = 'OK'

        def after_model(self, *args):
                pymol.cmd.do('delete bestmodel')
                chaincodes = pymol.cmd.get_names(enabled_only=1)
                #[os.rename(f, f.replace('input-alignment.ali.BL', 'model.BL')) for f in os.listdir('.') if not f.startswith('.')]


                try:
                        os.system('mv input.seg fasta.seg')
                except:
                        pass
                dope = self.get_best_model()
                pymol.cmd.do('load %s.pdb, bestmodel; as cartoon; alignto bestmodel'%milestones[3])
                pymol.cmd.do('load %s.pdb' %parameters[4])
                pymol.cmd.do('save bestmodel.pdb, bestmodel')
                # pymol.cmd.do('as cartoon')
                # pymol.cmd.do('alignto bestmodel')
                print "Generating energy profiles for templates... "
                chaincodes = pymol.cmd.get_names(enabled_only=1)
                for code in chaincodes:
                        os.system('python src/eval.py %s %s | tee profile.log' %(milestones[3],code))
                os.system('python src/profile_align.py %s %s | tee align.log' % (milestones[3], parameters[0]))
                milestones[5] = 'OK'
                print "Done!"

                # Get pdb parameters

                print " Extracting PDB info from bestmodel..."
                sleep(0.5)
                os.system('python src/tools/pdb_param.py %s.pdb | tee param.log'%milestones[3]) #get pdb parameters
                sleep(0.5)
                with open('param.log') as param:
                        mw, pI = param.readlines()[0].split() # add parameters to readable variables
                mw = str(round(float(mw)/1000, 2)) # convert MW to KDa format with 2 decimals
                pI = round(float(pI), 2)
                dope = round(float(dope), 2)
                residues = len(pymol.cmd.get_model("bestmodel and poly").get_residues()) # count residues using pymol API
                builder.get_object('entry9').set_text(str(residues))
                builder.get_object('entry10').set_text(str(mw)+' KDa')
                builder.get_object('entry11').set_text(str(pI))
                builder.get_object('entry12').set_text(str(dope))

                print 'Done!'




        def get_best_model(self, *args):

                total = ''
                data  = []
                with open('model.log','r') as log:
                        total = log.readlines()
                for line in total:
                        if 'target.B' in line:
                                data.append(line.split())
                        if 'non-standart' in line:
                                builder.get_object('label41').set_text("This is usually caused by non-standard residues, such \nas ligands, or by PDB files with missing atoms.")
                                builder.get_object('warning_dialog1').show()

                data_c = []
                for ent in data:
                        if 'target' in ent[0]:
                                data_c.append(ent)
                data_c.sort(key=lambda x: float(x[2]))

                print 'Adding best model to state file...'
                milestones[3]= data_c[0][0].split('.pdb')[0] # add the best model name to the milestones variable
                print '\nAll done!'
                return float(data_c[0][2])



        def plot_dope(self, *args):
                if '.B' not in milestones[3] or milestones[3]=='bestmodel':
                        self.get_best_model()

                from matplotlib.ticker import MaxNLocator

                print 'Ploting modeling DOPE graph'
                total = ''
                data  = []
                with open('model.log','r') as log:
                        total = log.readlines()
                for line in total:
                        if 'target.B' in line:
                                data.append(line.split())
                data_c = []
                for ent in data:
                        if 'target' in ent[0]:
                                data_c.append(ent)
                data_c
                vals = []

                for dope in data_c:
                        vals.append(float(dope[2]))

                win = Gtk.Window()
                win.connect("delete-event", Gtk.main_quit )
                win.set_default_size(600,450)
                win.set_title("Modeling DOPE graph")

                f = Figure(figsize=(5,4))
                a = f.add_subplot(1,1,1)
                x = range(1,len(vals)+1)


                xa = a.axes.get_xaxis()
                xa.set_major_locator(MaxNLocator(integer=True))
                a.set_ylabel('\nDOPE')
                a.set_xlabel('Residue\n')

                a.plot(x,vals, 'o-')

                vbox = Gtk.VBox()
                win.add(vbox)

                # Add canvas to vbox
                canvas = FigureCanvas(f)  # a Gtk.DrawingArea
                vbox.pack_start(canvas, True, True, 0)

                # Create toolbar
                toolbar = NavigationToolbar(canvas, win)
                vbox.pack_start(toolbar, False, False, 0)

                win.show_all()
                Gtk.main()


        def open_in_pymol(self, *args):
                if '.B' not in milestones[3] or milestones[3]=='bestmodel':
                        self.get_best_model()
                pymol.cmd.do('load %s.pdb, bestmodel'%milestones[3])

        def dope_per_residue(self, *args):
                if milestones[5] == 'evaluation':
                        os.system('python src/eval.py %s %s | tee profile.log' %(milestones[3],parameters[4]))
                        os.system('python src/profile_align.py %s %s | tee align.log' % (milestones[3], parameters[0]))
                        milestones[5] = 'OK'
	        e = modeller.environ()
	        al = modeller.alignment(e, file='input_alignment2d.ali')
	        code_model  = parameters[4]
	        code_model2 = parameters[5]
	        code_model3 = parameters[6]
	        code_templ  = milestones[3]
	        template    = self.get_profile('%s.profile' % code_model, al[code_model])
		rangelist   = template
                if parameters[5] != '':
			try:
			        template2   = self.get_profile('%s.profile' % code_model2, al[code_model2])
		                x3 = range(1, len(template2)+1)
				rangelist += template2
			except:
				pass


	        model       = self.get_profile('%s.profile' % code_templ, al['target'])
                rangelist += model
                rangelist = [x for x in rangelist if x is not None]

		mini = min(rangelist)
		maxi = max(rangelist)

                win = Gtk.Window()
                win.connect("delete-event", Gtk.main_quit )
                win.set_default_size(700,500)
                win.set_title("DOPE per residue for best model")


                f = plt.figure(figsize=(8,6), dpi=100)
                a = f.add_subplot(111)

                x1 = range(1, len(model)+1)
                x2 = range(1, len(template)+1)


                if x1>x2:
                        x = x1
                else:
                        x = x2


                a.plot(x1, model, color= 'red', linewidth=3, label = 'Best Model')
                a.plot(x2,template,  label = code_model)

                if parameters[5] != '':
			try:
				if x1>x3:
		                        x = x1
		                else:
		                        x = x3
	                        a.plot(x3,template2,  label = code_model2)
			except:
				pass

                if '/' in code_model3:

                        for codes in code_model3.split('/'):
                                print codes
			        try:

			                template3   = self.get_profile('%s.profile' % codes, al[codes])
                                        print 'o codes e', codes
		                        x3 = range(1, len(template3)+1)
				        rangelist += template3

		                        if x1>x3:
		                                x = x1
		                        else:
		                                x = x3
	                                a.plot(x3,template3,  label = codes)
			        except:
				        pass





                a.legend()
                a.set_xlim(1,len(x))
                a.set_ylim(1.15*float(mini),0.85*float(maxi))
                a.set_ylabel('DOPE')
                a.set_xlabel('Residue')
                a.set_title ('DOPE Profile per residue\n')

                vbox = Gtk.VBox()
                win.add(vbox)

                # Add canvas to vbox
                canvas = FigureCanvas(f)  # a Gtk.DrawingArea
                vbox.pack_start(canvas, True, True, 0)

                # Create toolbar
                toolbar = NavigationToolbar(canvas, win)
                vbox.pack_start(toolbar, False, False, 0)


                win.show_all()
                Gtk.main()



	def dope_per_residue_ref(self, *args):
		if '.BL00' not in milestones[6]:
			print "You don't have any refined model yet. Please run refinement first"
			builder.get_object('label41').set_text("You don't have any refined model yet. Please run refinement first")
			builder.get_object('warning_dialog1').show()

		else:

#	                os.system('python src/eval.py %s %s | tee profile_ref.log' %(milestones[3],milestones[6]))



			e = modeller.environ()
	        	a = modeller.alignment(e, file='profile.ali')
			code_model  = milestones[6]
			code_templ  = milestones[3]
			template    = self.get_profile('%s.profile' % code_model, a['target'])
			model       = self.get_profile('%s.profile' % code_templ, a['target'])


		        win = Gtk.Window()
		        win.connect("delete-event", Gtk.main_quit )
		        win.set_default_size(700,500)
		        win.set_title("DOPE per residue for best model")


		        f = plt.figure(figsize=(8,6), dpi=100)
		        a = f.add_subplot(111)

		        x1 = range(1, len(model)+1)
		        x2 = range(1, len(template)+1)


		        if x1>x2:
		                x = x1
		        else:
		                x = x2


		        a.plot(x1, model, color= 'red', label = 'Basic Model')
		        a.plot(x2,template, color= 'blue', linewidth=2,  label = 'Loop Ref Model')


		        a.legend()
		        a.set_xlim(1,len(x))
    			template = [x for x in template if x is not None]
    			a.set_ylim(1.15*float(min(template)),0.85*float(max(template)))
		        a.set_ylabel('DOPE')
		        a.set_xlabel('Residue')
		        a.set_title ('DOPE Profile per residue\n')

		        vbox = Gtk.VBox()
		        win.add(vbox)

		        # Add canvas to vbox
		        canvas = FigureCanvas(f)  # a Gtk.DrawingArea
		        vbox.pack_start(canvas, True, True, 0)

		        # Create toolbar
		        toolbar = NavigationToolbar(canvas, win)
		        vbox.pack_start(toolbar, False, False, 0)


		        win.show_all()
		        Gtk.main()




        def torsion(self, pdb_code):
                import math
                import Bio.PDB
                def degrees(rad_angle) :
                        if rad_angle is None :
                                return None
                        angle = rad_angle * 180 / math.pi
                        while angle > 180 :
                                angle = angle - 360
                        while angle < -180 :
                                angle = angle + 360
                        return angle

                structure = Bio.PDB.PDBParser().get_structure(pdb_code, "%s.pdb" % pdb_code)
                output    = []
                #Extract phi and psi angles
                for model in structure :
                        for chain in model :
                                polypeptides = Bio.PDB.CaPPBuilder().build_peptides(chain)
                                for poly_index, poly in enumerate(polypeptides) :
                                        phi_psi = poly.get_phi_psi_list()
                                        for res_index, residue in enumerate(poly) :
                                                phi, psi = phi_psi[res_index]
                                                if phi and psi :
                                                        #Don't write output when missing an angle

                                                        output.append(("%s:Chain%s:%s%i\t%f\t%f\n" \
                                                                % (pdb_code, str(chain.id), residue.resname,
                                                                   residue.id[1], degrees(phi), degrees(psi),
                                                                   )).split('\t'))


                phi_scat = []
                psi_scat = []
                resname  = []

                for i in output:
                        phi_scat.append(i[1])
                        psi_scat.append(i[2])
                        resname.append(i[0].split(':')[-1])

                #convert strings to floats for ploting function
                map(float, phi_scat)
                map(float, psi_scat)

                return resname, phi_scat, psi_scat


        def ramachandran_chooser(self, *args):
                builder.get_object('dialog1').show()


        def ramachandran_close(self, *args):
                builder.get_object('dialog1').hide()


        def ramachandran(self, *args):
                global rama_pdb


                if builder.get_object('radiobutton3').get_active():
                        pdb = milestones[3]
                        if '.B' not in milestones[3] or milestones[3]=='bestmodel':
                                self.get_best_model()

                        if milestones[3] not in pymol.cmd.get_names('all'):
                                pymol.cmd.do('load %s.pdb, bestmodel' %milestones[3])

                        if milestones[3] == 'bestmodel':
                                builder.get_object('label41').set_text('No best model found.\nTry to import a modeling session if you\nalready modeled your protein')
                                builder.get_object('warning_dialog1').show()
                        rama_pdb = 'bestmodel'

                elif builder.get_object('radiobutton4').get_active():
                        pdb = milestones[6]
                        if '.B' not in milestones[6] or milestones[6]=='ref_bestmodel':
                                self.get_best_loop_ref()

                        if milestones[6] not in pymol.cmd.get_names('all'):
                                pymol.cmd.do('load %s.pdb, ref_bestmodel' %milestones[6])
                                pymol.cmd.do('as cartoon')
                        if milestones[6] == 'bestmodel':
                                builder.get_object('label41').set_text('No best model found.\nTry to import a modeling session if you\nalready modeled your protein')
                                builder.get_object('warning_dialog1').show()
                        rama_pdb = 'ref_bestmodel'

                if '.pdb' in pdb:
                        pdb = pdb.split('.pdb')[0]

                #Retrieve torsion angles from best model
                resname, phi, psi = self.torsion(pdb)
                builder.get_object('dialog1').hide()
                self.plot_picker(phi,psi,resname)


        def plot_picker(self, phi, psi, resname):
            import numpy as np
            import re
            if rama_pdb == 'bestmodel':
                    pdb = milestones[3]
            elif rama_pdb == 'ref_bestmodel':
                    pdb = milestones[6]
            elif rama_pdb == 'selectedPDB_ramachandran':
                    pdb = builder.get_object('filechooserbutton2').get_filename().split('.pdb')[0]
            rama_sel = rama_pdb
            print 'Running Procheck....'
            try:
                    os.system('$HOME/Vitamin/src/procheck/procheck.scr %s.pdb 2' %pdb)
            except:
                    os.system('procheck %s.pdb 2' %pdb)

            if milestones[4] != rama_pdb:

                    milestones[4] = rama_pdb

            try:    #create folder for results
                    os.stat('procheck')
                    #os.system('rm procheck/*')
            except:
                    os.mkdir('procheck')

            savefiles = ['*.sdh', '*.sco', 'tplot.log', '*.new', 'anglen.log', 'procheck.prm', '*.out', '*.pln', '*.lan', '*.nb', 'pplot.log', 'bplot.log', 'secstr.log', '*.ps', '*.rin', '*.sum', 'clean.log', 'nb.log']

            print 'Moving files to procheck folder...'


            for files in savefiles:         #move every file create since last session to results
                    os.system('mv %s procheck/' %files)
            print "Done."

            summary = open('procheck/%s.sum'%pdb).read()
            print summary

            sumlist = open('procheck/%s.sum'%pdb).readlines()

            for line in sumlist:
                    if "Ramachandran plot:" in line:
                            rama_stat = line.split()


            class PointBrowser:
                    """
                    Click on a point to select and highlight it -- the data that
                    generated the point will be shown in the lower axes.  Use the 'n'
                    and 'p' keys to browse through the next and previous points
                    """
                    def __init__(self):
                        self.lastind = 0

                        self.selected,  = ax.plot([xs[0]], [ys[0]], 'o', ms=12, alpha=0.4,
                                                  color='yellow', visible=False)

                    def onpress(self, event):
                        if self.lastind is None: return
                        if event.key not in ('n', 'p'): return
                        if event.key=='n': inc = 1
                        else:  inc = -1


                        self.lastind += inc
                        self.lastind = np.clip(self.lastind, 0, len(xs)-1)
                        self.update()

                    def onpick(self, event):

                       if event.artist!=line: return True

                       N = len(event.ind)
                       if not N: return True

                       # the click locations
                       x = event.mouseevent.xdata
                       y = event.mouseevent.ydata


                       distances = np.hypot(x-float(xs[event.ind]), y-float(ys[event.ind]))
                       indmin = distances.argmin()
                       dataind = event.ind[indmin]

                       self.lastind = dataind
                       self.update()

                    def update(self):
                        if self.lastind is None: return

                        dataind = self.lastind


                        #self.selected.set_visible(True)
                        self.selected.set_data(xs[dataind], ys[dataind])

                        #self.text.set_text('selected: %d'%dataind)

                        ano = plt.annotate(
                                resname[dataind],
                                xy = (xs[dataind], ys[dataind]), xytext = (-20, 20),
                                textcoords = 'offset points', ha = 'right', va = 'bottom',
                                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
                        #show and center residue in pymol
                        pymol.cmd.do('select resi %s' %resindex[dataind])
                        pymol.cmd.do('label (n. CA and resi %s), resn+resi' %resindex[dataind])
                        pymol.cmd.do('show sticks, resi %s and %s' %(resindex[dataind],rama_sel))
                        pymol.cmd.do('color red, resi %s and %s' %(resindex[dataind],rama_sel))
                        pymol.cmd.do('center resi %s and %s' %(resindex[dataind],rama_sel))
                        fig.canvas.draw()


            if __name__ == '__main__':

                        #separate resn and resi
                        r = re.compile("([a-zA-Z]+)([0-9]+)")
                        resindex = [int(r.match(res).groups()[1]) for res in resname]

                        win = Gtk.Window()
                        win.connect("delete-event", Gtk.main_quit )
                        win.set_default_size(950,900)
                        win.set_title("Ramachandran Plots")

                        xs = phi
                        ys = psi

                        # define subplot and aspect ratio
                        fig, (ax) = plt.subplots(1)
                        ax.set(aspect=1)

                        #ax.plot (phi, psi,  'g.')
                        line, = ax.plot(xs, ys, 'o', picker=5)  # 5 points tolerance



                        #define plot area properties
                        ax.axis ([-180,180,-180,180])
                        ax.imshow (plt.imread('gui/regions.png')   , aspect='equal', extent = (-180,180,-180,180))
                        ax.set_xlabel(r'$\phi$', fontsize=40)
                        ax.set_ylabel(r'$\psi$', fontsize=40)
                        ax.set_title ('General\n', fontsize=20)
                        ax.set_xticks(np.arange(-180,180+1, 30))
                        ax.set_yticks(np.arange(-180,180+1, 30))


                        browser = PointBrowser()

                        vbox = Gtk.VBox()

                        win.add(vbox)
                        # Create a Gtk.DrawingArea
                        canvas = FigureCanvas(fig)
                        vbox.pack_start(canvas, True, True, 0)
                        # Create label
                        label = Gtk.Label()
                        label.set_markup('<big>Favoured: %s    Allowed: %s    Generously Allowed: %s    Disallowed: %s</big>' %(rama_stat[3],rama_stat[5],rama_stat[7],rama_stat[9]))
                        vbox.pack_start(label, False, False,0)

                        # Create toolbar
                        toolbar = NavigationToolbar(canvas, win)
                        vbox.pack_start(toolbar, False, False, 0)


                        fig.canvas.mpl_connect('pick_event', browser.onpick)
                        fig.canvas.mpl_connect('key_press_event', browser.onpress)




                        win.show_all()
                        Gtk.main()






        def plot_picker_standalone(self,*args):
            import numpy as np
            import re

            filechooser = builder.get_object('filechooserdialog2')
            filechooser.hide()
            pdb = filechooser.get_filenames()[0]
            pymol.cmd.do('load %s'%pdb)
            pymol.cmd.do('as cartoon')
            pymol.cmd.do('color green')


            if '.pdb' in pdb:
                        pdb = pdb.split('.pdb')[0]

            resname, phi, psi = self.torsion(pdb)

            print 'Running Procheck....'
            try:
                    os.system('$HOME/Vitamin/src/procheck/procheck.scr %s.pdb 2' %pdb)
            except:
                    os.system('procheck %s.pdb 2' %pdb)

            pdb=pdb.split('/')[-1]

            try:    #create folder for results
                    os.stat('procheck_standalone')

            except:
                    os.mkdir('procheck_standalone')

            savefiles = ['*.sdh', '*.sco', 'tplot.log', '*.new', 'anglen.log', 'procheck.prm', '*.out', '*.pln', '*.lan', '*.nb', 'pplot.log', 'bplot.log', 'secstr.log', '*.ps', '*.rin', '*.sum', 'clean.log', 'nb.log']

            print 'Moving files to procheck folder...'


            for files in savefiles:         #move every file create since last session to results
                    os.system('mv %s procheck_standalone/' %files)
            print "Done."

            summary = open('procheck_standalone/%s.sum'%pdb).read()
            print summary

            sumlist = open('procheck_standalone/%s.sum'%pdb).readlines()

            for line in sumlist:
                    if "Ramachandran plot:" in line:
                            rama_stat = line.split()


            class PointBrowser:
                    """
                    Click on a point to select and highlight it -- the data that
                    generated the point will be shown in the lower axes.  Use the 'n'
                    and 'p' keys to browse through the next and previous points
                    """
                    def __init__(self):
                        self.lastind = 0

                        self.selected,  = ax.plot([xs[0]], [ys[0]], 'o', ms=12, alpha=0.4,
                                                  color='yellow', visible=False)

                    def onpress(self, event):
                        if self.lastind is None: return
                        if event.key not in ('n', 'p'): return
                        if event.key=='n': inc = 1
                        else:  inc = -1


                        self.lastind += inc
                        self.lastind = np.clip(self.lastind, 0, len(xs)-1)
                        self.update()

                    def onpick(self, event):

                       if event.artist!=line: return True

                       N = len(event.ind)
                       if not N: return True

                       # the click locations
                       x = event.mouseevent.xdata
                       y = event.mouseevent.ydata


                       distances = np.hypot(x-float(xs[event.ind]), y-float(ys[event.ind]))
                       indmin = distances.argmin()
                       dataind = event.ind[indmin]

                       self.lastind = dataind
                       self.update()

                    def update(self):
                        if self.lastind is None: return

                        dataind = self.lastind


                        #self.selected.set_visible(True)
                        self.selected.set_data(xs[dataind], ys[dataind])

                        #self.text.set_text('selected: %d'%dataind)

                        ano = plt.annotate(
                                resname[dataind],
                                xy = (xs[dataind], ys[dataind]), xytext = (-20, 20),
                                textcoords = 'offset points', ha = 'right', va = 'bottom',
                                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
                        #show and center residue in pymol
                        pymol.cmd.do('select resi %s' %resindex[dataind])
                        pymol.cmd.do('label (n. CA and resi %s), resn+resi' %resindex[dataind])
                        pymol.cmd.do('show sticks, resi %s' %(resindex[dataind]))
                        pymol.cmd.do('color red, resi %s' %(resindex[dataind]))
                        pymol.cmd.do('center resi %s' %(resindex[dataind]))
                        fig.canvas.draw()


            if __name__ == '__main__':

                        #separate resn and resi
                        r = re.compile("([a-zA-Z]+)([0-9]+)")
                        resindex = [int(r.match(res).groups()[1]) for res in resname]

                        win = Gtk.Window()
                        win.connect("delete-event", Gtk.main_quit )
                        win.set_default_size(950,900)
                        win.set_title("Ramachandran Plots")

                        xs = phi
                        ys = psi

                        # define subplot and aspect ratio
                        fig, (ax) = plt.subplots(1)
                        ax.set(aspect=1)

                        #ax.plot (phi, psi,  'g.')
                        line, = ax.plot(xs, ys, 'o', picker=5)  # 5 points tolerance



                        #define plot area properties
                        ax.axis ([-180,180,-180,180])
                        ax.imshow (plt.imread('gui/regions.png')   , aspect='equal', extent = (-180,180,-180,180))
                        ax.set_xlabel(r'$\phi$', fontsize=40)
                        ax.set_ylabel(r'$\psi$', fontsize=40)
                        ax.set_title ('General\n', fontsize=20)
                        ax.set_xticks(np.arange(-180,180+1, 30))
                        ax.set_yticks(np.arange(-180,180+1, 30))


                        browser = PointBrowser()

                        vbox = Gtk.VBox()

                        win.add(vbox)
                        # Create a Gtk.DrawingArea
                        canvas = FigureCanvas(fig)
                        vbox.pack_start(canvas, True, True, 0)
                        # Create label
                        label = Gtk.Label()
                        label.set_markup('<big>Favoured: %s    Allowed: %s    Generously Allowed: %s    Disallowed: %s</big>' %(rama_stat[3],rama_stat[5],rama_stat[7],rama_stat[9]))
                        vbox.pack_start(label, False, False,0)

                        # Create toolbar
                        toolbar = NavigationToolbar(canvas, win)
                        vbox.pack_start(toolbar, False, False, 0)


                        fig.canvas.mpl_connect('pick_event', browser.onpick)
                        fig.canvas.mpl_connect('key_press_event', browser.onpress)




                        win.show_all()
                        Gtk.main()




        def loop_selector(self, *args):

                if 'bestmodel' not in pymol.cmd.get_names('all'):
                        self.open_in_pymol()
                pymol.cmd.do('show cartoon')
#
#               Check for loaded best model
                if milestones[5] == 'evaluation':
                        os.system('python src/eval.py %s %s | tee profile.log' %(milestones[3],parameters[4]))
                        os.system('python src/profile_align.py %s %s | tee align.log' % (milestones[3], parameters[0])) #create alignment with itself for profile plot

                        milestones[5] = 'OK'
	        e = modeller.environ()
	        a = modeller.alignment(e, file='profile.ali')
	        code = milestones[3]
	        model = self.get_profile(milestones[3]+'.profile', a[code])



#               Plot graph with span select and show in pymol
                import numpy as np
                from matplotlib.widgets import SpanSelector
                def onselect(xmin, xmax):
                    indmin, indmax = np.searchsorted(x, (xmin, xmax))
                    indmax = min(len(x)-1, indmax)

                    thisx = x[indmin:indmax]


                    f.canvas.draw()
                    builder.get_object('entry13').set_text(str(round(xmin))[:-2])
                    builder.get_object('entry14').set_text(str(round(xmax))[:-2])

                    self.zoom_selected_loop()


                win = Gtk.Window()
                win.connect("delete-event", Gtk.main_quit )
                win.set_default_size(700,550)
                win.set_title("Select Loop Region")


                f = plt.figure(figsize=(8,6), dpi=100)
                a = f.add_subplot(111)

                x = range(1, len(model)+1)

                a.plot(x, model, color= 'red',  label = 'Best Model')
		a.axhline(y=-0.03, linewidth=1, color='black')
                a.legend()
                a.set_xlim(1,len(x))
                a.set_ylabel('DOPE')
                a.set_xlabel('Residue')
                a.set_title ('DOPE Profile per residue\n')

                vbox = Gtk.VBox()
                win.add(vbox)

                # Add canvas to vbox
                canvas = FigureCanvas(f)  # a Gtk.DrawingArea
                vbox.pack_start(canvas, True, True, 0)

                # Create toolbar
                toolbar = NavigationToolbar(canvas, win)
                vbox.pack_start(toolbar, False, False, 0)




                # set useblit True on gtkagg for enhanced performance
                span = SpanSelector(a, onselect, 'horizontal', useblit=False,
                                    rectprops=dict(alpha=0.5, facecolor='red') )


                win.show_all()
                Gtk.main()


        def zoom_selected_loop(self, *args):
                res_ini = builder.get_object('entry13').get_text()
                res_fin = builder.get_object('entry14').get_text()

                pymol.cmd.do('as cartoon')
                pymol.cmd.do('color green')
                if res_ini != '' and res_fin != '':
                        pymol.cmd.do('select resi %s-%s' % (res_ini, res_fin))

                pymol.cmd.do('color red, sele')
                pymol.cmd.do('orient sele; delete sele')


        def loop_ref(self, *args):
                ini = builder.get_object('entry13').get_text()
                fin = builder.get_object('entry14').get_text()
                num = builder.get_object('entry15').get_text()
                models_pymol = filter( lambda s: ('ref_bestmodel' in s), pymol.cmd.get_names(enabled_only=1))

                # models_pymol = [s for s in models_pymol if 'ref_' in s]
                if not len(models_pymol):
                    models_pymol = ['bestmodel']
                import re
                def natural_sort(l):
                    convert = lambda text: int(text) if text.isdigit() else text.lower()
                    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
                    return sorted(l, key = alphanum_key)

                pymol.cmd.do('save %s.pdb, %s' %(models_pymol[-1],models_pymol[-1]))

                self.start_spinner('Refining models')

                def thread_run(*args):
                        model_ret = os.system('python src/loop_model.py %s %s %s %s %s | tee loop_ref.log' % (ini, fin, models_pymol[-1], num, parameters[7].split('-')[0]))
                        self.get_ref_profiles()
                        GObject.idle_add(cleanup_func, model_ret)
                def cleanup_func(model_ret):
                        Gdk.threads_enter()
                        self.spinnerdialog.destroy()
                        Gdk.threads_leave()
                        self.after_ref_model()


                t = threading.Thread(target=thread_run)
                t.start()


        def after_ref_model(self,*args):
                #pymol.cmd.do('delete ref_bestmodel')
                models_pymol = pymol.cmd.get_names(enabled_only=1)
                dope = self.get_best_loop_ref()
                pymol.cmd.do('load %s.pdb, ref_bestmodel%s; as cartoon; alignto bestmodel'%(milestones[6],' '.join(models_pymol).count('ref_')))
                try:
                        os.system('python src/eval.py %s %s | tee profile_ref.log' %(milestones[3],milestones[6]))
                except:
                        print "Your bestmodel or ref_bestmodel PDB is missing from your folder."

                # pymol.cmd.do('as cartoon')
                # pymol.cmd.do('alignto bestmodel')
                print "Generating energy profiles... "
                sleep(0.5)
                print "Done!"



        def get_best_loop_ref(self, *args):

                #Compute modeling dope profile
                Dope = imp.load_source('get_dope','src/dope_ref.py')
                results = Dope.get_dope()
                codes = []

                for dope in results:
                        if len(dope)>2:
                                codes.append(dope.split()[0])
                        else:
                                pass


                vals = []

                for dope in results:
                        try:
                                vals.append(dope.split()[1])
                        except:
                                pass

                milestones[6] = codes[vals.index(min(vals))].split('.pdb')[0]
                print '\nAll done!'
                return max(vals)



#plotar o grafico dos molpdfs para os loops




        def plot_dope_loop_ref(self, *args):


                from matplotlib.ticker import MaxNLocator

                print 'Ploting modeling molpdf graph'
                Dope = imp.load_source('get_dope','src/dope_ref.py')
                results = Dope.get_dope()
                vals = []

                for dope in results:
                        try:
                                vals.append(dope.split()[1])
                        except:
                                pass

                win = Gtk.Window()
                win.connect("delete-event", Gtk.main_quit )
                win.set_default_size(600,450)
                win.set_title("Refinement MOLPDF graph")

                f = Figure(figsize=(5,4))
                a = f.add_subplot(1,1,1)
                x = range(1,len(vals)+1)


                xa = a.axes.get_xaxis()
                xa.set_major_locator(MaxNLocator(integer=True))
                a.set_ylabel('MOLPDF')
                a.set_xlabel('Residue')
                a.plot(x,vals, 'o-')

                vbox = Gtk.VBox()
                win.add(vbox)

                # Add canvas to vbox
                canvas = FigureCanvas(f)  # a Gtk.DrawingArea
                vbox.pack_start(canvas, True, True, 0)

                # Create toolbar
                toolbar = NavigationToolbar(canvas, win)
                vbox.pack_start(toolbar, False, False, 0)

                win.show_all()
                Gtk.main()


        def get_ref_profiles(self, *args):

                ref_pdbs = glob.glob("*.BL*.pdb")
                for pdb in ref_pdbs:
                        os.system('python src/eval_ref.py %s' % pdb)


        def show_refloop_pymol(self, *args):
                print 'Loading structures in PyMOL...'
                ref_pdbs = glob.glob("*.BL*.pdb")
                for pdb in ref_pdbs:
                        pymol.cmd.load(pdb)
                sleep(1)
                try:
                        pymol.cmd.do('alignto %s' % ref_pdbs[0])
                        pymol.cmd.do('as cartoon')
                except:
                        pass


        def preview_ref_changes(self, *args):
                pymol.cmd.do('alignto %s' % ref_pdbs[0])
                try:
                        pymol.cmd.do('select resi %s-%s' % (res_ini, res_fin) )
                except:
                        print "Please set loop residues from DOPE Graph"

        def plot_dope_ref_model(self,*args):


	        e = modeller.environ()
	        a = modeller.alignment(e, file='profile.ali')

                win = Gtk.Window()
                win.connect("delete-event", Gtk.main_quit )
                win.set_default_size(700,550)
                win.set_title("DOPE profile of refined loops")


                f = plt.figure(figsize=(8,6), dpi=100)
                ax = f.add_subplot(111)
                model = self.get_profile(milestones[3]+'.profile', a['target'])
                x = range(1, len(model)+1)
                ax.plot(x, model , linewidth=1,label = 'Best Model')

                for pdb in glob.glob("*.BL*.pdb"):
        	        model = self.get_profile(pdb+'.profile', a['target'])
                        ax.plot(x, model,  label = pdb)

                ax.legend()
                ax.axhline(y=-0.03, linewidth=1, color='black')

                ax.set_xlim(1,len(x))
                ax.set_ylabel('DOPE')
                ax.set_xlabel('Residue')
                ax.set_title ('DOPE Profile per residue\n')

                vbox = Gtk.VBox()
                win.add(vbox)

                # Add canvas to vbox
                canvas = FigureCanvas(f)  # a Gtk.DrawingArea
                vbox.pack_start(canvas, True, True, 0)

                # Create toolbar
                toolbar = NavigationToolbar(canvas, win)
                vbox.pack_start(toolbar, False, False, 0)




                win.show_all()
                Gtk.main()





#########################################################
#########################################################
#########################################################
#########################################################

        def create_new_folder(*args):
            pass






import imp
import sys, os

# autocompletion
import readline
import rlcompleter
readline.parse_and_bind('tab: complete')

# pymol environment
moddir='/usr/share/pymol/'
sys.path.insert(0, moddir)
os.environ['PYMOL_PATH'] = os.path.join(moddir, '')

# pymol launching
import pymol
pymol.finish_launching()


def vitamin_spash(*args):
    os.system('clear')
    print ""+'\n'
    print "   VITAMIN - Visual Table Modeller Interface, version 1.0.0"
    print "    Created by Gustavo M. A. de Lima and Fernando V. Maluf"
    print ""
    print "    ViTaMIn is open-source software under LGPL 3 licence."
    print "         Please consider giving helpful feedback"
    print ""




os.system('clear')


builder = Gtk.Builder()
builder.add_from_file("vit.glade")
builder.get_object('main').show_all()
builder.connect_signals(Handler())
vitamin_spash()






class ThreadClass(threading.Thread):
	def run(self):
               	Gtk.main()

ThreadClass().start()

#
#
# import Tkinter
#
# root = Tkinter.Tk()
#
# def autosave():
#     os.system('rm config')
#     while True:
#         # do something you want
#         open('config','w').write(','.join(parameters)+'\n'+','.join(milestones))
#         time.sleep(30)
#
# saver = threading.Thread(target=autosave)
# saver.start()
