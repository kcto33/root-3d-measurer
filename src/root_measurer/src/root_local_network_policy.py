"""Native measurement computes from local files; deny sockets in its process."""
import sys

_installed=False
DENIED=('socket.connect','socket.sendto','socket.getaddrinfo','socket.gethostbyname','socket.gethostbyaddr')


def disable_measurement_network():
    global _installed
    if not _installed:
        def guard(event,args):
            if event in DENIED:raise RuntimeError('NATIVE_MEASUREMENT_OFFLINE_NETWORK_DISABLED')
        sys.addaudithook(guard)
        _installed=True
    return network_policy_state()


def network_policy_state():
    return dict(active=_installed,policy='deny measurement socket connections and DNS',
                socket_events_denied=list(DENIED) if _installed else [],system_network_changed=False)
