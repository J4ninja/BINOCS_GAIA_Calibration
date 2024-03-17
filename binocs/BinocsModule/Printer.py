
class Printer():

    def print_initial_results(self, options, mag, info, summary):
        print("    Writing intial results to '%s--binary.txt'" % (options['data']))
        out = open(options['data']+"--binary.txt", "w")
        for s in range(mag.shape[0]):
            # Print out star to file
            outstr = "%16s %9.5f %9.5f " % (info[s][0], info[s][1], info[s][2])
            for i in range(17): outstr += "%6.3f " % mag[s,2*i]
            outstr += "   %2d %2d %6.3f %6.4f %6.3f %6.4f %5.2f   %6.3f %6.4f %5.2f" % (summary[s,8], 0, summary[s,0], summary[s,1], summary[s,2], summary[s,3], summary[s,4], summary[s,5], summary[s,6], summary[s,7])
            print(outstr, file=out)
        out.close()

    def print_synthetic_results(self, options, synth, binary, synth_summary):
        print("    Writing synthetic results to '%s--synth.txt'" % (options['data']))
        out = open(options['data']+"--synth.txt", "w")
        for s in range(synth.shape[0]):
            # Print out star to file
            outstr = "%9.5f %9.5f " % (binary[s,0], binary[s,1])
            for i in range(17): outstr += "%6.3f " % synth[s,2*i]
            outstr += "   %6.3f %6.4f %6.3f %6.4f %5.2f   %6.3f %6.4f %5.2f" % (synth_summary[s,0], synth_summary[s,1], synth_summary[s,2], synth_summary[s,3], synth_summary[s,4], synth_summary[s,5], synth_summary[s,6], synth_summary[s,7])
            print(outstr, file=out)
        out.close()

    def print_minimum_mass_ratios(self, options, minq_synth, minq_dm, minq_mod):
        print("    Writing minimum mass ratio results to '%s--minq.txt'" % (options['data']))
        out = open(options['data']+"--minq.txt", "w")
        for m in range(len(minq_synth)):
            print("%5.2f  %5.3f  %5.3f  %5.3f" % (m*minq_dm, minq_synth[m], minq_mod[m], max([minq_synth[m], minq_mod[m], 0.3])), file=out)
        out.close()

    def print_final_results(self, options, mag, summary, minmass, minq_synth, minq_dm, info):
        print("    Writing updated results to '%s--results.txt'" % (options['data']))
        out = open(options['data']+"--results.txt", 'w')
        for s in range(mag.shape[0]):
            # Determine new binary flag and masses
            if summary[s,0] == 0:
                # This star is a non-member
                mpri, mprie, msec, msece = 0, 0, 0, 0
                bflag = summary[s,8]
            elif summary[s,2] < minmass or summary[s,2] / summary[s,0] < minq_synth[int(summary[s,0]//minq_dm)] or summary[s,8] == 1:
                # This star is a single
                mpri, mprie, msec, msece = summary[s,5], summary[s,6], 0, 0
                bflag = 1
            else:
                # This star is a binary
                mpri, mprie, msec, msece = summary[s,0], summary[s,1], summary[s,2], summary[s,3]
                bflag = 2
            # Print out star to file
            outstr = "%16s %9.5f %9.5f " % (info[s][0], info[s][1], info[s][2])
            for i in range(17): outstr += "%6.3f " % mag[s,2*i]
            outstr += "   %2d %2d %6.3f %6.4f %6.3f %6.4f" % (bflag, info[s][3], mpri, mprie, msec, msece)
            print(outstr, file=out)
        out.close()	