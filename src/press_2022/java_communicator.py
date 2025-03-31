from py4j.java_gateway import (
    JavaGateway, CallbackServerParameters, GatewayParameters,
    launch_gateway)
from py4j.java_collections import ListConverter

from generate_barcodes import pickle_barcodes, load_barcodes
from decode_reads import move_to_gpu, get_reads, call_reads
from calc_triage_recall import get_index_recall, get_triage_positions
from util import *

import struct  # TODO


class JavaCommunicator(object):

    def __init__(self):
        self.barcodes = None
        self.barcode_dict_cpu = None
        self.barcode_dict_gpu = None
        self.valid_barcodes = None
        self.gateway = None

    def setBarcodes(self, barcodes_dataset):
        # Slow because of low throughput of p4j
        self.barcodes = barcodes_dataset

        if barcodes_dataset is not None:
            codes = [encode(code.toString()[:]) for code in barcodes_dataset.getSequences()]
            self.barcode_dict_cpu = pickle_barcodes(codes)
            self.barcode_dict_gpu = move_to_gpu(self.barcode_dict_cpu)

    def setBarcodesFile(self, barcodes_filename):
        self.barcode_dict_cpu = load_barcodes(barcodes_filename)
        self.barcode_dict_gpu = move_to_gpu(self.barcode_dict_cpu)

    @measure_execution_time
    def setValidBarcodes(self, valid_barcodes):
        self.valid_barcodes = valid_barcodes

        if valid_barcodes is not None:
            barcode_dict_cpu_filtered = filter_barcode_dict(self.barcode_dict_cpu, valid_barcodes)
            self.barcode_dict_gpu = move_to_gpu(barcode_dict_cpu_filtered)
        else:
            self.barcode_dict_gpu = move_to_gpu(self.barcode_dict_cpu)

    def calcCosineIndex(self, read_seq, k, valid_barcodes=None):
        if valid_barcodes is not None:
            self.setValidBarcodes(valid_barcodes)

        if isinstance(read_seq, str):
            seq = read_seq
        else:
            seq = array(encode(read_seq.toString()[:]))

        N = self.barcode_dict_gpu["N"]
        M = self.barcode_dict_gpu["M"]

        dists = torch.zeros(N, dtype=torch.float, device=device)  # will hold distances for each read
        allrank = torch.zeros(N, dtype=torch.float, device=device)

        # primary and secondary triage
        mer = seqtomer(seq)

        cosvec = torch.zeros(64, dtype=torch.float, device=device)
        cosvec[mer] = self.barcode_dict_gpu["coses"][k-1, torch.arange(len(mer), dtype=torch.long, device=device)]

        dists[:] = torch.sum(torch.unsqueeze(cosvec, dim=0) * self.barcode_dict_gpu["cosvecs"][k-1],
                             dim=1)  # all cosine distances

        allrank[:] = torch.argsort(dists[:], descending=True)
        # Convert torch array to python list
        if self.valid_barcodes is not None:
            candidates = [self.valid_barcodes[int(i)] for i in allrank]
        else:
            candidates = [int(i) for i in allrank]
        # Convert to Java list
        java_list = ListConverter().convert(candidates, gateway._gateway_client)

        return java_list

    def calcHammingIndex(self, read_seq, valid_barcodes=None):
        if valid_barcodes is not None:
            self.setValidBarcodes(valid_barcodes)

        if isinstance(read_seq, str):
            seq = read_seq
        else:
            seq = array(encode(read_seq.toString()[:]))

        N = self.barcode_dict_gpu["N"]
        M = self.barcode_dict_gpu["M"]

        dists = torch.zeros(N, dtype=torch.float, device=device)  # will hold distances for each read
        allrank = torch.zeros(N, dtype=torch.float, device=device)

        # primary and secondary triage
        mer = seqtomer(seq)
        foo = int64(uint64(mertobin(mer)))  # need to cast 64 bits to a type known to torch
        seqbin = torch.tensor(foo, dtype=torch.int64, device=device)
        xored = torch.bitwise_xor(seqbin, self.barcode_dict_gpu["allbitmaps"])
        dists[:] = mypopcount(xored)  # all Hamming distances

        allrank[:] = torch.argsort(dists[:], descending=False)
        # Convert torch array to python list
        if self.valid_barcodes is not None:
            candidates = [self.valid_barcodes[int(i)] for i in allrank]
        else:
            candidates = [int(i) for i in allrank]
        # Convert to Java list
        java_list = ListConverter().convert(candidates, gateway._gateway_client)

        return java_list

    @measure_execution_time
    def getRecallAtPos(self, reads_filename, labels_filename):
        reads = get_reads(reads_filename)
        answers = fromfile(labels_filename, dtype=int32)

        recall = get_index_recall(self.barcode_dict_gpu, reads, answers)
        java_array = gateway.new_array(gateway.jvm.double, recall.shape[0])
        for i, rec in enumerate(recall):
            java_array[i] = double(rec)

        return java_array

    @measure_execution_time
    def writeFrequencyAtPos(self, reads_filename, labels_filename, out_filename):
        reads = get_reads(reads_filename)
        answers = fromfile(labels_filename, dtype=int32)
        positions = get_triage_positions(self.barcode_dict_gpu, reads, answers)
        positions.tofile(out_filename)

    @measure_execution_time
    def getTriageCandidates(self, reads_filename):
        torch.set_grad_enabled(False)
        reads = get_reads(reads_filename)

        N = self.barcode_dict_gpu["N"]
        M = self.barcode_dict_gpu["M"]
        r = reads.shape[0]

        Ncos = self.barcode_dict_gpu["cosvecs"].shape[0]
        dists = torch.zeros(Ncos + 1, N, dtype=torch.float, device=device)  # will hold distances for each read
        allrank = torch.zeros(Ncos + 1, N, dtype=torch.int32, device=device)

        candidates_lists = zeros(((Ncos + 1) * r, N), dtype=int32)

        for j, seq in enumerate(reads[:r]):
            # primary and secondary triage
            mer = seqtomer(seq)
            foo = int64(uint64(mertobin(mer)))  # need to cast 64 bits to a type known to torch
            seqbin = torch.tensor(foo, dtype=torch.int64, device=device)
            xored = torch.bitwise_xor(seqbin, self.barcode_dict_gpu["allbitmaps"])
            dists[0, :] = 64. - mypopcount(xored)  # all Hamming distances
            cosvec = torch.zeros(Ncos, 64, dtype=torch.float, device=device)
            for k in range(Ncos):
                cosvec[k, mer] = self.barcode_dict_gpu["coses"][
                    k, torch.arange(len(mer), dtype=torch.long, device=device)]
            dists[1:, :] = torch.sum(torch.unsqueeze(cosvec, dim=1) * self.barcode_dict_gpu["cosvecs"],
                                     dim=2)  # all cosine distances
            for k in range(Ncos + 1):
                allrank[k, :] = torch.argsort(dists[k, :], descending=True)  # rank them all
                candidates_lists[j*(Ncos+1) + k] = allrank[k].cpu().numpy()

        return candidates_lists

    @measure_execution_time
    def writeTriageCandidates(self, reads_filename, out_filename, valid_barcodes_filename=None):
        if valid_barcodes_filename is not None:
            valid_barcodes = fromfile(valid_barcodes_filename, dtype=int32)
            self.setValidBarcodes(valid_barcodes)

        candidates_lists = self.getTriageCandidates(reads_filename)
        # Write arrays to binary file (ints are 4 byte with byte[0] + byte[1]*2^8 + byte[3]*2^16 ...)
        candidates_lists.tofile(out_filename)

    @measure_execution_time
    def decodeReads(self, reads_filename, decodes_filename, exec_times_filename=None, dists_filename=None, pos=100, L=8):
        reads = get_reads(reads_filename)
        torch.set_grad_enabled(False)

        print(exec_times_filename)  # TODO

        best, best_dists, exec_time_triages, exec_time_levenshtein = call_reads(reads, self.barcode_dict_gpu, pos=pos, L=L)

        print(exec_time_triages, exec_time_levenshtein) # TODO

        best.tofile(decodes_filename)
        if dists_filename is not None:
            best_dists.tofile(dists_filename)

        if exec_times_filename is not None:
            with open(exec_times_filename, "w") as writer:
                writer.write(f"{exec_time_triages}\n{exec_time_levenshtein}")
                print(f"{exec_time_triages}\n{exec_time_levenshtein}")  # TODO

    class Java:
        implements = ["de.uni_halle.barcode_calling.matchers.indexes.PressCommunicator"]


# Global methods
def filter_barcode_dict(barcode_dict, valid_barcodes):
    Nnew = len(valid_barcodes)
    M = barcode_dict["M"]

    allseqs = []
    alltrimers = []
    allbitmaps = zeros(Nnew, dtype=uint64)
    cosvecs = torch.zeros((N_COS, Nnew, 64), dtype=torch.float)
    coses = barcode_dict["coses"]

    for i, barcode_id in enumerate(valid_barcodes):
        allseqs.append(barcode_dict["allseqs"][barcode_id])
        alltrimers.append(barcode_dict["alltrimers"][barcode_id])
        allbitmaps[i] = barcode_dict["allbitmaps"][barcode_id]

        for k in range(N_COS):
            cosvecs[k, i] = barcode_dict["cosvecs"][k][barcode_id]

    pickledict = {"N": Nnew, "M": M, "allseqs": allseqs,
                  "alltrimers": alltrimers, "allbitmaps": allbitmaps, "coses": coses, "cosvecs": cosvecs}

    return pickledict


# Main
if __name__ == "__main__":
    print("Listener started and listens!")

    communicator = JavaCommunicator()

    gateway = JavaGateway(
        gateway_parameters=GatewayParameters(auto_convert=True),
        python_server_entry_point=communicator)

    gateway.start_callback_server()

    """
    port = launch_gateway()

    # connect python side to Java side with Java dynamic port and start python
    # callback server with a dynamic port
    gateway = JavaGateway(
        gateway_parameters=GatewayParameters(port=port),
        callback_server_parameters=CallbackServerParameters(port=0),
        python_server_entry_point=communicator)

    # retrieve the port on which the python callback server was bound to.
    python_port = gateway.get_callback_server().get_listening_port()

    # tell the Java side to connect to the python callback server with the new
    # python port. Note that we use the java_gateway_server attribute that
    # retrieves the GatewayServer instance.
    gateway.java_gateway_server.resetCallbackClient(
        gateway.java_gateway_server.getCallbackClient().getAddress(),
        python_port)
    """

